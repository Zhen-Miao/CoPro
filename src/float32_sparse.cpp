#include <Rcpp.h>

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>
#include <thread>
#include <unordered_map>
#include <vector>

using namespace Rcpp;

namespace {

struct CellKeyF32 {
  std::int64_t x;
  std::int64_t y;
  std::int64_t z;

  bool operator==(const CellKeyF32& other) const {
    return x == other.x && y == other.y && z == other.z;
  }
};

struct CellKeyF32Hash {
  std::size_t operator()(const CellKeyF32& key) const {
    std::size_t seed = std::hash<std::int64_t>()(key.x);
    seed ^= std::hash<std::int64_t>()(key.y) + 0x9e3779b9U +
      (seed << 6) + (seed >> 2);
    seed ^= std::hash<std::int64_t>()(key.z) + 0x9e3779b9U +
      (seed << 6) + (seed >> 2);
    return seed;
  }
};

struct FloatEdge {
  int column;
  float distance;
};

std::int64_t checked_grid_size_f32(double span, double radius) {
  const double value = std::floor(span / radius) + 1.0;
  const double max_index = static_cast<double>(
    std::numeric_limits<std::int64_t>::max()
  );
  if (!R_finite(value) || value < 1.0 || value > max_index) {
    stop("Grid extent is too large for the fixed-radius neighbor search.");
  }
  return static_cast<std::int64_t>(value);
}

CellKeyF32 point_cell_f32(
    const NumericMatrix& coordinates,
    int row,
    const std::vector<double>& origin,
    const std::vector<std::int64_t>& grid_size,
    double radius) {
  std::int64_t cell[3] = {0, 0, 0};
  const int dimensions = coordinates.ncol();
  for (int axis = 0; axis < dimensions; ++axis) {
    double raw = std::floor(
      (coordinates(row, axis) - origin[axis]) / radius
    );
    if (!R_finite(raw)) {
      stop("Coordinates must contain only finite values.");
    }
    raw = std::max(0.0, raw);
    raw = std::min(raw, static_cast<double>(grid_size[axis] - 1));
    cell[axis] = static_cast<std::int64_t>(raw);
  }
  return CellKeyF32{cell[0], cell[1], cell[2]};
}

bool inside_f32(std::int64_t value, std::int64_t size) {
  return value >= 0 && value < size;
}

// The packed payload is stored little-endian. On the little-endian platforms R
// targets (x86-64, ARM64, Windows) that is native byte order, so a single
// aligned load/store is both correct and reads/writes the value in one machine
// instruction -- the same contiguous-typed-array access that makes Matrix's
// double path fast, instead of assembling each value from four separate bytes.
// R allocates vector data on at least an 8-byte boundary, so the 4-byte access
// is always aligned. The explicit byte path is retained for a big-endian build
// so the on-buffer format stays little-endian and remains portable.
#if defined(_WIN32) || defined(__LITTLE_ENDIAN__) || \
    (defined(__BYTE_ORDER__) && __BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__)
#  define COPRO_F32_NATIVE_LE 1
#endif

inline void write_float(Rbyte* destination, R_xlen_t index, float value) {
#ifdef COPRO_F32_NATIVE_LE
  std::memcpy(destination + index * 4, &value, sizeof(value));
#else
  std::uint32_t bits;
  std::memcpy(&bits, &value, sizeof(bits));
  const R_xlen_t offset = index * 4;
  destination[offset] = bits & 0xffu;
  destination[offset + 1] = (bits >> 8) & 0xffu;
  destination[offset + 2] = (bits >> 16) & 0xffu;
  destination[offset + 3] = (bits >> 24) & 0xffu;
#endif
}

inline float read_float(const Rbyte* source, R_xlen_t index) {
#ifdef COPRO_F32_NATIVE_LE
  float value;
  std::memcpy(&value, source + index * 4, sizeof(value));
  return value;
#else
  const R_xlen_t offset = index * 4;
  const std::uint32_t bits =
    static_cast<std::uint32_t>(source[offset]) |
    (static_cast<std::uint32_t>(source[offset + 1]) << 8) |
    (static_cast<std::uint32_t>(source[offset + 2]) << 16) |
    (static_cast<std::uint32_t>(source[offset + 3]) << 24);
  float value;
  std::memcpy(&value, &bits, sizeof(value));
  return value;
#endif
}

float type7_quantile(
    std::vector<float>& values,
    double probability,
    int repetitions = 1) {
  if (values.empty()) stop("Cannot compute a quantile of no values.");
  if (repetitions < 1) stop("Quantile repetitions must be positive.");
  if (values.size() == 1U) return values[0];

  const double represented_size =
    static_cast<double>(values.size()) * repetitions;
  const double h = (represented_size - 1.0) *
    probability;
  const std::size_t lower = static_cast<std::size_t>(std::floor(h));
  const double fraction = h - static_cast<double>(lower);
  const std::size_t lower_source =
    lower / static_cast<std::size_t>(repetitions);
  std::nth_element(
    values.begin(), values.begin() + lower_source, values.end()
  );
  const float lower_value = values[lower_source];
  if (fraction == 0.0 ||
      lower + 1U >= static_cast<std::size_t>(represented_size)) {
    return lower_value;
  }
  const std::size_t upper_source =
    (lower + 1U) / static_cast<std::size_t>(repetitions);
  if (upper_source == lower_source) return lower_value;
  std::nth_element(
    values.begin(), values.begin() + upper_source, values.end()
  );
  const float upper_value = values[upper_source];
  return static_cast<float>(
    static_cast<double>(lower_value) +
      fraction * (
        static_cast<double>(upper_value) -
          static_cast<double>(lower_value)
      )
  );
}

// Grid buckets in compressed form.
//
// The obvious structure is unordered_map<CellKeyF32, vector<int>>, but that is
// one heap-allocated vector per occupied grid cell -- on the order of one per
// point at the radii CoPro uses -- each carrying 24 bytes of book-keeping plus
// its own allocation, and scattering the point lists all over the heap. Here
// the map holds a single int per cell and the membership lists live in one
// contiguous counting-sorted array, so iterating a bucket is a linear walk.
//
// Members stay in increasing point order, matching the push_back order of the
// structure this replaces, so enumeration order (and therefore output) is
// unchanged.
struct FlatGridBuckets {
  std::unordered_map<CellKeyF32, int, CellKeyF32Hash> slot;
  std::vector<int> start;
  std::vector<int> items;

  bool empty_cell(const CellKeyF32& key, int* begin, int* end) const {
    const auto found = slot.find(key);
    if (found == slot.end()) return true;
    *begin = start[found->second];
    *end = start[found->second + 1];
    return false;
  }
};

FlatGridBuckets build_flat_buckets(
    const NumericMatrix& B,
    int n_b,
    const std::vector<double>& origin,
    const std::vector<std::int64_t>& grid_size,
    double radius) {
  FlatGridBuckets grid;
  grid.slot.reserve(static_cast<std::size_t>(n_b * 1.3) + 1U);
  std::vector<int> cell_of(static_cast<std::size_t>(n_b));
  std::vector<int> counts;

  for (int row = 0; row < n_b; ++row) {
    const CellKeyF32 key = point_cell_f32(B, row, origin, grid_size, radius);
    const auto inserted =
      grid.slot.insert(std::make_pair(key, static_cast<int>(counts.size())));
    if (inserted.second) counts.push_back(0);
    const int id = inserted.first->second;
    cell_of[row] = id;
    ++counts[id];
  }

  grid.start.assign(counts.size() + 1U, 0);
  for (std::size_t id = 0; id < counts.size(); ++id) {
    grid.start[id + 1] = grid.start[id] + counts[id];
  }
  grid.items.resize(static_cast<std::size_t>(n_b));
  std::vector<int> cursor(grid.start.begin(), grid.start.end() - 1);
  for (int row = 0; row < n_b; ++row) {
    grid.items[cursor[cell_of[row]]++] = row;
  }
  return grid;
}

template <typename Callback>
std::size_t visit_neighbors(
    int row_begin,
    int row_end,
    const NumericMatrix& A,
    const NumericMatrix& B,
    const std::vector<CellKeyF32>& cells_a,
    const FlatGridBuckets& buckets,
    const std::vector<std::int64_t>& grid_size,
    double radius_squared,
    double percentile,
    double scaling_factor,
    bool truncate_low_distance,
    bool symmetric,
    Callback callback) {
  const int dimensions = A.ncol();
  std::size_t count = 0U;
  for (int a = row_begin; a < row_end; ++a) {
    const CellKeyF32 base = cells_a[a];
    for (int dz = dimensions == 3 ? -1 : 0;
         dz <= (dimensions == 3 ? 1 : 0); ++dz) {
      const std::int64_t z = base.z + dz;
      if (dimensions == 3 && !inside_f32(z, grid_size[2])) continue;
      for (int dy = -1; dy <= 1; ++dy) {
        const std::int64_t y = base.y + dy;
        if (!inside_f32(y, grid_size[1])) continue;
        for (int dx = -1; dx <= 1; ++dx) {
          const std::int64_t x = base.x + dx;
          if (!inside_f32(x, grid_size[0])) continue;
          int bucket_begin = 0;
          int bucket_end = 0;
          if (buckets.empty_cell(CellKeyF32{x, y, z},
                                 &bucket_begin, &bucket_end)) {
            continue;
          }
          for (int slot = bucket_begin; slot < bucket_end; ++slot) {
            const int b = buckets.items[slot];
            // A symmetric within-type kernel stores its strict upper
            // triangle. Self-pairs are excluded, matching .frnnGrid(A, NULL).
            if (symmetric && b <= a) continue;
            double distance_squared = 0.0;
            for (int axis = 0; axis < dimensions; ++axis) {
              const double delta = A(a, axis) - B(b, axis);
              distance_squared += delta * delta;
            }
            if (distance_squared <= radius_squared) {
              double distance = std::sqrt(distance_squared);
              if (truncate_low_distance && distance < percentile) {
                distance = percentile;
              }
              callback(
                a, b,
                static_cast<float>(distance * scaling_factor)
              );
              ++count;
            }
          }
        }
      }
    }
  }
  return count;
}

int normalized_thread_count(int requested, int tasks) {
  if (requested < 1) stop("n_threads must be at least one.");
  return std::max(1, std::min(requested, tasks));
}

// Repack an R column-major double matrix as row-major float32.
//
// The sparse operators walk a kernel row's nonzeros and, for each one, read all
// n_columns_out values of one row of X. In R's column-major layout those values
// sit n_rows doubles apart, so a single nonzero touches one cache line per
// column -- for the nPCA = 10-40 matrices CoPro uses, that is 10-40 lines (and
// often as many pages) per nonzero. Paying one O(n * ncol) pass up front to
// make each X row contiguous turns the inner loop into a unit-stride sweep.
//
// The float32 conversion is the same static_cast the inner loop performed on
// each access, just hoisted, so every product is formed from an identical
// value and the accumulation order is untouched: results are bit-identical.
// The buffer costs n * ncol * 4 bytes (24 MB at 200k x 30) and is shared
// read-only across workers.
std::vector<float> pack_row_major_f32(
    const double* source, int n_rows_in, int n_columns_in) {
  std::vector<float> packed(
    static_cast<std::size_t>(n_rows_in) *
      static_cast<std::size_t>(n_columns_in)
  );
  for (int column = 0; column < n_columns_in; ++column) {
    const double* input_column =
      source + static_cast<std::size_t>(n_rows_in) * column;
    for (int row = 0; row < n_rows_in; ++row) {
      packed[
        static_cast<std::size_t>(row) *
          static_cast<std::size_t>(n_columns_in) + column
      ] = static_cast<float>(input_column[row]);
    }
  }
  return packed;
}

void validate_csr(
    const IntegerVector& p,
    const IntegerVector& j,
    const RawVector& x,
    const IntegerVector& dims,
    bool symmetric = false) {
  if (dims.size() != 2) stop("dims must have length two.");
  if (dims[0] < 0 || dims[1] < 0) stop("dims must be nonnegative.");
  if (p.size() != static_cast<R_xlen_t>(dims[0]) + 1) {
    stop("CSR row pointer length does not match nrow.");
  }
  if (j.size() > std::numeric_limits<int>::max()) {
    stop("CSR nonzero count exceeds the 32-bit index limit.");
  }
  if (p[0] != 0 || p[p.size() - 1] != j.size()) {
    stop("Invalid CSR row pointer.");
  }
  for (R_xlen_t index = 0; index < p.size(); ++index) {
    if (p[index] < 0 || p[index] > j.size() ||
        (index > 0 && p[index] < p[index - 1])) {
      stop("CSR row pointers must be nonnegative and non-decreasing.");
    }
  }
  for (R_xlen_t index = 0; index < j.size(); ++index) {
    if (j[index] < 0 || j[index] >= dims[1]) {
      stop("CSR column index is out of bounds.");
    }
  }
  if (x.size() != static_cast<R_xlen_t>(j.size()) * 4) {
    stop("Float32 value buffer does not match the CSR nonzero count.");
  }
  if (symmetric && dims[0] != dims[1]) {
    stop("A symmetric CSR kernel must be square.");
  }
  if (symmetric) {
    for (int row = 0; row < dims[0]; ++row) {
      for (int position = p[row]; position < p[row + 1]; ++position) {
        if (j[position] <= row) {
          stop("A symmetric CSR kernel must store only its strict upper triangle.");
        }
      }
    }
  }
}

}  // namespace


// Build all requested Gaussian kernels for one cross-cell-type block.
// Distances are calculated in double and immediately retained as float32 in a
// row-compressed temporary. Persistent kernel values are emitted directly as
// IEEE-754 float32 bytes; no double sparse matrix is materialized.
//
// [[Rcpp::export(rng = false)]]
List float32_csr_gaussian_kernels_cpp(
    NumericMatrix A,
    NumericMatrix B,
    NumericVector sigmas,
    double percentile,
    double scaling_factor,
    double lower_limit,
    double upper_quantile,
    bool truncate_low_distance = true,
    bool symmetric = false,
    int normalization = 0,
    int n_threads = 1) {
  const int n_a = A.nrow();
  const int n_b = B.nrow();
  const int dimensions = A.ncol();
  if (dimensions != 2 && dimensions != 3) {
    stop("Float32 fixed-radius kernels support only 2-D or 3-D coordinates.");
  }
  if (B.ncol() != dimensions) {
    stop("A and B must have the same coordinate dimension.");
  }
  if (symmetric && n_a != n_b) {
    stop("Symmetric float32 kernels require equally sized coordinate blocks.");
  }
  if (n_a == 0 || n_b == 0) stop("Coordinate blocks must be nonempty.");
  if (sigmas.size() == 0) stop("At least one sigma is required.");
  if (!R_finite(percentile) || percentile < 0.0) {
    stop("percentile must be nonnegative and finite.");
  }
  if (!R_finite(scaling_factor) || scaling_factor <= 0.0) {
    stop("scaling_factor must be positive and finite.");
  }
  if (!R_finite(lower_limit) ||
      lower_limit <= 0.0 || lower_limit >= 1.0) {
    stop("lower_limit must lie strictly between zero and one.");
  }
  if (!R_finite(upper_quantile) ||
      upper_quantile <= 0.0 || upper_quantile >= 1.0) {
    stop("upper_quantile must lie strictly between zero and one.");
  }
  if (normalization < 0 || normalization > 3) {
    stop("normalization must be 0 (none), 1 (global), 2 (row), or 3 (column).");
  }

  double max_sigma = 0.0;
  for (R_xlen_t index = 0; index < sigmas.size(); ++index) {
    if (!R_finite(sigmas[index]) || sigmas[index] <= 0.0) {
      stop("All sigma values must be positive and finite.");
    }
    max_sigma = std::max(max_sigma, sigmas[index]);
  }
  const double support_multiplier = std::sqrt(-2.0 * std::log(lower_limit));
  const double radius = support_multiplier * max_sigma / scaling_factor *
    (1.0 + 1e-6);
  const double radius_squared = radius * radius;

  std::vector<double> origin(dimensions);
  std::vector<std::int64_t> grid_size(dimensions);
  for (int axis = 0; axis < dimensions; ++axis) {
    double lower = A(0, axis);
    double upper = A(0, axis);
    if (!R_finite(lower)) stop("Coordinates must be finite.");
    for (int row = 1; row < n_a; ++row) {
      const double value = A(row, axis);
      if (!R_finite(value)) stop("Coordinates must be finite.");
      lower = std::min(lower, value);
      upper = std::max(upper, value);
    }
    for (int row = 0; row < n_b; ++row) {
      const double value = B(row, axis);
      if (!R_finite(value)) stop("Coordinates must be finite.");
      lower = std::min(lower, value);
      upper = std::max(upper, value);
    }
    origin[axis] = lower;
    grid_size[axis] = checked_grid_size_f32(upper - lower, radius);
  }

  std::vector<CellKeyF32> cells_a;
  cells_a.reserve(n_a);
  for (int row = 0; row < n_a; ++row) {
    cells_a.push_back(
      point_cell_f32(A, row, origin, grid_size, radius)
    );
  }
  const FlatGridBuckets buckets =
    build_flat_buckets(B, n_b, origin, grid_size, radius);

  // Estimate the edge count from a row sample to avoid vector-capacity
  // doubling at large n while retaining a single enumeration pass.
  const int sample_rows = std::min(n_a, 2048);
  const std::size_t sample_count = visit_neighbors(
    0, sample_rows, A, B, cells_a, buckets, grid_size,
    radius_squared, percentile, scaling_factor,
    truncate_low_distance, symmetric,
    [](int, int, float) {}
  );
  const double estimated = sample_rows > 0 ?
    static_cast<double>(sample_count) *
      static_cast<double>(n_a) / static_cast<double>(sample_rows) :
    0.0;
  const double reserve_double = std::max(
    static_cast<double>(n_a) * 8.0,
    std::ceil(estimated * 1.15)
  );
  if (reserve_double > static_cast<double>(std::numeric_limits<int>::max())) {
    stop("A single float32 sparse kernel block would exceed the 32-bit compressed-index limit.");
  }

  // Enumerate neighbours in parallel over disjoint row ranges. Each worker
  // fills its own edge buffer and its own slice of the row pointer; the buffers
  // are then concatenated in thread order. Because the ranges are contiguous
  // and ascending, and each worker walks its rows in order, the concatenation
  // reproduces exactly the row-major edge order the serial loop produced --
  // the CSR output is unchanged, not merely equivalent.
  //
  // No R API call may occur off the main thread, so checkUserInterrupt() cannot
  // run inside the workers. It is checked before and after instead; the phase
  // is correspondingly shorter, and the per-sigma loop below still polls.
  const int enumeration_threads =
    normalized_thread_count(n_threads, std::max(1, n_a));
  std::vector<std::vector<FloatEdge> > thread_edges(enumeration_threads);
  std::vector<std::size_t> edge_pointer(
    static_cast<std::size_t>(n_a) + 1U, 0U
  );
  std::vector<std::size_t> thread_zero_counts(enumeration_threads, 0U);
  std::vector<float> thread_minimum(
    enumeration_threads, std::numeric_limits<float>::infinity());
  std::vector<int> range_begin(enumeration_threads + 1);
  for (int index = 0; index <= enumeration_threads; ++index) {
    range_begin[index] = static_cast<int>(
      static_cast<std::int64_t>(n_a) * index / enumeration_threads);
  }

  checkUserInterrupt();
  {
    std::vector<std::thread> workers;
    workers.reserve(enumeration_threads);
    std::atomic<bool> enumeration_failed(false);
    for (int index = 0; index < enumeration_threads; ++index) {
      workers.emplace_back([&, index]() {
        try {
          const int row_begin = range_begin[index];
          const int row_end = range_begin[index + 1];
          std::vector<FloatEdge>& local = thread_edges[index];
          local.reserve(static_cast<std::size_t>(
            reserve_double * (row_end - row_begin) / std::max(1, n_a)) + 16U);
          std::size_t zeros = 0U;
          float minimum = std::numeric_limits<float>::infinity();
          for (int row = row_begin; row < row_end; ++row) {
            visit_neighbors(
              row, row + 1, A, B, cells_a, buckets, grid_size,
              radius_squared, percentile, scaling_factor,
              truncate_low_distance, symmetric,
              [&](int, int column, float distance) {
                if (distance == 0.0f) {
                  ++zeros;
                } else {
                  minimum = std::min(minimum, distance);
                }
                local.push_back(FloatEdge{column, distance});
              }
            );
            // Per-row counts now; turned into global offsets after the join.
            edge_pointer[static_cast<std::size_t>(row) + 1U] =
              local.size();
          }
          thread_zero_counts[index] = zeros;
          thread_minimum[index] = minimum;
        } catch (...) {
          enumeration_failed.store(true);
        }
      });
    }
    for (std::thread& worker : workers) worker.join();
    if (enumeration_failed.load()) {
      stop("A float32 neighbour-enumeration worker failed.");
    }
  }
  checkUserInterrupt();

  std::size_t zero_distance_count = 0U;
  float minimum_nonzero = std::numeric_limits<float>::infinity();
  std::size_t total_edges = 0U;
  std::size_t private_capacity_bytes = 0U;
  for (int index = 0; index < enumeration_threads; ++index) {
    zero_distance_count += thread_zero_counts[index];
    minimum_nonzero = std::min(minimum_nonzero, thread_minimum[index]);
    total_edges += thread_edges[index].size();
    private_capacity_bytes +=
      thread_edges[index].capacity() * sizeof(FloatEdge);
  }
  if (total_edges > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
    stop("A single float32 sparse kernel block exceeds the 32-bit compressed-index limit.");
  }

  // Concatenating k private buffers into one contiguous array costs a transient
  // second copy: `edges` is reserved at full size while the buffers that feed
  // it are still live. Buffers are released as they are consumed, so the excess
  // decays over the loop, but the peak is reached at the first insert and is
  // real -- roughly 8 bytes per edge on top of the array itself. With one
  // thread there is nothing to concatenate, so move the single buffer instead
  // and pay nothing; that keeps the memory-tightest configuration at the same
  // peak as the serial code this replaced. See RESULTS.md for the measured
  // multi-thread figure.
  std::vector<FloatEdge> edges;
  std::size_t peak_edge_bytes = private_capacity_bytes;
  if (enumeration_threads == 1) {
    edges = std::move(thread_edges[0]);
    // A single range already carries offsets relative to the whole array.
  } else {
    edges.reserve(total_edges);
    // The private buffers and the array they feed are all live at the first
    // insert. That instant is the high-water mark of the whole build, and it
    // is what `temporary_bytes` reports -- the capacity retained afterwards is
    // strictly smaller and would understate the parallel path's cost.
    peak_edge_bytes =
      private_capacity_bytes + total_edges * sizeof(FloatEdge);
    std::size_t offset = 0U;
    for (int index = 0; index < enumeration_threads; ++index) {
      std::vector<FloatEdge>& local = thread_edges[index];
      edges.insert(edges.end(), local.begin(), local.end());
      // Each worker recorded counts relative to its own buffer; shift them onto
      // the concatenated one.
      for (int row = range_begin[index]; row < range_begin[index + 1]; ++row) {
        edge_pointer[static_cast<std::size_t>(row) + 1U] +=
          offset;
      }
      offset += local.size();
      std::vector<FloatEdge>().swap(local);  // release as we go
    }
  }

  if (zero_distance_count > 0U) {
    warning("Zero distances detected in a float32 kernel block; applying dense-compatible nearest-distance handling.");
    for (FloatEdge& edge : edges) {
      if (edge.distance == 0.0f) {
        // If no nonzero pair lies inside the retained support, the true
        // block-global minimum lies outside it and would be dropped by every
        // requested kernel. Infinity represents that drop without an O(n*m)
        // fallback scan.
        edge.distance = std::isfinite(static_cast<double>(minimum_nonzero)) ?
          minimum_nonzero : std::numeric_limits<float>::infinity();
      }
    }
  }

  List kernels(sigmas.size());
  NumericVector nonzeros(sigmas.size());
  NumericVector stored_nonzeros(sigmas.size());
  for (R_xlen_t sigma_index = 0;
       sigma_index < sigmas.size(); ++sigma_index) {
    checkUserInterrupt();
    const float inverse_sigma = static_cast<float>(1.0 / sigmas[sigma_index]);
    std::vector<float> quantile_values;
    quantile_values.reserve(edges.size());
    IntegerVector pointer(n_a + 1);

    for (int row = 0; row < n_a; ++row) {
      int row_count = 0;
      for (std::size_t position = edge_pointer[row];
           position < edge_pointer[row + 1]; ++position) {
        const float scaled = edges[position].distance * inverse_sigma;
        const float weight = std::exp(-0.5f * scaled * scaled);
        if (weight >= static_cast<float>(lower_limit)) {
          quantile_values.push_back(weight);
          ++row_count;
        }
      }
      if (pointer[row] > std::numeric_limits<int>::max() - row_count) {
        stop("A float32 sparse kernel exceeds the 32-bit index limit.");
      }
      pointer[row + 1] = pointer[row] + row_count;
    }

    const int raw_nnz = pointer[n_a];
    nonzeros[sigma_index] = raw_nnz;
    const bool expand_symmetric =
      symmetric && (normalization == 2 || normalization == 3);
    const bool output_symmetric = symmetric && !expand_symmetric;
    if (raw_nnz == 0) {
      stored_nonzeros[sigma_index] = 0;
      kernels[sigma_index] = List::create(
        _["p"] = pointer,
        _["j"] = IntegerVector(0),
        _["x"] = RawVector(0),
        _["Dim"] = IntegerVector::create(n_a, n_b),
        _["transposed"] = false,
        _["symmetric"] = output_symmetric
      );
      continue;
    }

    const float upper_clip = type7_quantile(
      quantile_values, upper_quantile, symmetric ? 2 : 1
    );
    quantile_values.clear();
    quantile_values.shrink_to_fit();

    // Normalization sums and scaling factors deliberately remain float32.
    // Frobenius diagnostics use a separate double-accumulation path.
    std::vector<float> row_scale(
      normalization == 1 || normalization == 2 ?
        static_cast<std::size_t>(n_a) : 0U,
      0.0f
    );
    std::vector<float> column_scale(
      normalization == 3 ? static_cast<std::size_t>(n_b) : 0U,
      0.0f
    );
    if (normalization != 0) {
      for (int row = 0; row < n_a; ++row) {
        for (std::size_t position = edge_pointer[row];
             position < edge_pointer[row + 1]; ++position) {
          const float scaled = edges[position].distance * inverse_sigma;
          float weight = std::exp(-0.5f * scaled * scaled);
          if (weight < static_cast<float>(lower_limit)) continue;
          weight = std::min(weight, upper_clip);
          const int column = edges[position].column;
          if (normalization == 1 || normalization == 2) {
            row_scale[row] += weight;
            if (symmetric) row_scale[column] += weight;
          } else {
            column_scale[column] += weight;
            if (symmetric) column_scale[row] += weight;
          }
        }
      }

      if (normalization == 1) {
        std::vector<float> positive_row_sums;
        positive_row_sums.reserve(row_scale.size());
        for (const float value : row_scale) {
          if (value > 1e-5f) positive_row_sums.push_back(value);
        }
        float global_scale = 1.0f;
        if (!positive_row_sums.empty()) {
          const float median = type7_quantile(
            positive_row_sums, 0.5
          );
          if (std::isfinite(median) && median > 0.0f) {
            global_scale = 1.0f / median;
          }
        }
        std::fill(row_scale.begin(), row_scale.end(), global_scale);
      } else if (normalization == 2) {
        for (float& value : row_scale) {
          value = value > 1e-4f ? 1.0f / value : 1.0f;
        }
      } else {
        for (float& value : column_scale) {
          value = value > 1e-4f ? 1.0f / value : 1.0f;
        }
      }
    }

    const auto normalize_weight = [&](
        float weight, int row, int column) {
      if (normalization == 1 || normalization == 2) {
        return weight * row_scale[row];
      }
      if (normalization == 3) {
        return weight * column_scale[column];
      }
      return weight;
    };

    IntegerVector output_pointer;
    int output_nnz = raw_nnz;
    if (normalization == 0) {
      output_pointer = pointer;
    } else {
      std::vector<std::int64_t> row_counts(
        static_cast<std::size_t>(n_a), 0
      );
      for (int row = 0; row < n_a; ++row) {
        for (std::size_t position = edge_pointer[row];
             position < edge_pointer[row + 1]; ++position) {
          const float scaled = edges[position].distance * inverse_sigma;
          float weight = std::exp(-0.5f * scaled * scaled);
          if (weight < static_cast<float>(lower_limit)) continue;
          weight = std::min(weight, upper_clip);
          const int column = edges[position].column;
          if (normalize_weight(weight, row, column) >=
              static_cast<float>(lower_limit)) {
            ++row_counts[row];
          }
          if (expand_symmetric &&
              normalize_weight(weight, column, row) >=
                static_cast<float>(lower_limit)) {
            ++row_counts[column];
          }
        }
      }
      output_pointer = IntegerVector(n_a + 1);
      for (int row = 0; row < n_a; ++row) {
        const std::int64_t next =
          static_cast<std::int64_t>(output_pointer[row]) +
            row_counts[row];
        if (next > std::numeric_limits<int>::max()) {
          stop("A normalized float32 sparse kernel exceeds the 32-bit index limit.");
        }
        output_pointer[row + 1] = static_cast<int>(next);
      }
      output_nnz = output_pointer[n_a];
    }

    stored_nonzeros[sigma_index] = output_nnz;
    IntegerVector column_index(output_nnz);
    RawVector values(static_cast<R_xlen_t>(output_nnz) * 4);
    Rbyte* value_pointer = RAW(values);
    if (expand_symmetric) {
      std::vector<int> cursor(static_cast<std::size_t>(n_a));
      for (int row = 0; row < n_a; ++row) {
        cursor[row] = output_pointer[row];
      }
      for (int row = 0; row < n_a; ++row) {
        for (std::size_t position = edge_pointer[row];
             position < edge_pointer[row + 1]; ++position) {
          const float scaled = edges[position].distance * inverse_sigma;
          float weight = std::exp(-0.5f * scaled * scaled);
          if (weight < static_cast<float>(lower_limit)) continue;
          weight = std::min(weight, upper_clip);
          const int column = edges[position].column;
          const float forward = normalize_weight(weight, row, column);
          if (forward >= static_cast<float>(lower_limit)) {
            const int output_position = cursor[row]++;
            column_index[output_position] = column;
            write_float(value_pointer, output_position, forward);
          }
          const float reverse = normalize_weight(weight, column, row);
          if (reverse >= static_cast<float>(lower_limit)) {
            const int output_position = cursor[column]++;
            column_index[output_position] = row;
            write_float(value_pointer, output_position, reverse);
          }
        }
      }
    } else {
      int output_position = 0;
      for (int row = 0; row < n_a; ++row) {
        for (std::size_t position = edge_pointer[row];
             position < edge_pointer[row + 1]; ++position) {
          const float scaled = edges[position].distance * inverse_sigma;
          float weight = std::exp(-0.5f * scaled * scaled);
          if (weight < static_cast<float>(lower_limit)) continue;
          weight = std::min(weight, upper_clip);
          const int column = edges[position].column;
          weight = normalize_weight(weight, row, column);
          if (weight < static_cast<float>(lower_limit)) continue;
          column_index[output_position] = column;
          write_float(value_pointer, output_position, weight);
          ++output_position;
        }
      }
    }

    kernels[sigma_index] = List::create(
      _["p"] = output_pointer,
      _["j"] = column_index,
      _["x"] = values,
      _["Dim"] = IntegerVector::create(n_a, n_b),
      _["transposed"] = false,
      _["symmetric"] = output_symmetric
    );
  }

  return List::create(
    _["kernels"] = kernels,
    _["nonzeros"] = nonzeros,
    _["stored_nonzeros"] = stored_nonzeros,
    _["candidate_pairs"] = static_cast<double>(edges.size()),
    // Worst-case bound on edge storage allocated at one time, not the capacity
    // still held at return: with more than one enumeration thread the private
    // buffers and the concatenated array are briefly allocated together. It is
    // a bound rather than a measurement -- reserved pages are only made
    // resident as they are written, and the buffers are freed as they are
    // consumed, so measured RSS stays far below it (see RESULTS.md).
    _["temporary_bytes"] =
      static_cast<double>(peak_edge_bytes) +
      static_cast<double>(edge_pointer.capacity() * sizeof(std::size_t)),
    _["zero_distances_replaced"] =
      static_cast<double>(zero_distance_count)
  );
}


// Parallel float32 Y = X_left' K X_right for a row-compressed float32 K.
// Rows are disjoint across workers, so no atomics or n_rows x n_PC temporary
// are required. Thread-local PC x PC matrices are reduced after joining.
//
// [[Rcpp::export(rng = false)]]
NumericMatrix float32_csr_xky_cpp(
    IntegerVector p,
    IntegerVector j,
    RawVector x,
    IntegerVector dims,
    NumericMatrix x_left,
    NumericMatrix x_right,
    int n_threads = 1,
    bool symmetric = false) {
  validate_csr(p, j, x, dims, symmetric);
  const int n_rows = dims[0];
  const int n_columns = dims[1];
  if (x_left.nrow() != n_rows) {
    stop("x_left row count does not match the kernel.");
  }
  if (x_right.nrow() != n_columns) {
    stop("x_right row count does not match the kernel.");
  }
  const int p_left = x_left.ncol();
  const int p_right = x_right.ncol();
  const int threads = normalized_thread_count(n_threads, std::max(1, n_rows));

  const int* row_pointer = INTEGER(p);
  const int* column_index = INTEGER(j);
  const Rbyte* values = RAW(x);

  // Row-major float32 copies so each nonzero reads one contiguous run per
  // operand instead of p_left/p_right separate cache lines. See
  // pack_row_major_f32().
  const std::vector<float> left_packed =
    pack_row_major_f32(REAL(x_left), n_rows, p_left);
  const std::vector<float> right_packed =
    pack_row_major_f32(REAL(x_right), n_columns, p_right);
  const float* left_data = left_packed.data();
  const float* right_data = right_packed.data();

  std::vector<std::vector<float> > partial(
    threads,
    std::vector<float>(
      static_cast<std::size_t>(p_left) *
        static_cast<std::size_t>(p_right),
      0.0f
    )
  );
  std::vector<std::thread> workers;
  workers.reserve(threads);
  std::atomic<bool> failed(false);

  for (int thread_index = 0; thread_index < threads; ++thread_index) {
    const int row_begin = static_cast<int>(
      static_cast<std::int64_t>(n_rows) * thread_index / threads
    );
    const int row_end = static_cast<int>(
      static_cast<std::int64_t>(n_rows) * (thread_index + 1) / threads
    );
    workers.emplace_back([&, thread_index, row_begin, row_end]() {
      try {
        std::vector<float> kernel_row_right(p_right, 0.0f);
        std::vector<float> kernel_row_left(
          symmetric ? p_left : 0, 0.0f
        );
        std::vector<float>& output = partial[thread_index];
        for (int row = row_begin; row < row_end; ++row) {
          std::fill(
            kernel_row_right.begin(), kernel_row_right.end(), 0.0f
          );
          if (symmetric) {
            std::fill(
              kernel_row_left.begin(), kernel_row_left.end(), 0.0f
            );
          }
          for (int position = row_pointer[row];
               position < row_pointer[row + 1]; ++position) {
            const int column = column_index[position];
            const float weight = read_float(values, position);
            const float* right_row = right_data +
              static_cast<std::size_t>(column) *
              static_cast<std::size_t>(p_right);
            for (int b = 0; b < p_right; ++b) {
              kernel_row_right[b] += weight * right_row[b];
            }
            if (symmetric) {
              const float* left_row = left_data +
                static_cast<std::size_t>(column) *
                static_cast<std::size_t>(p_left);
              for (int a = 0; a < p_left; ++a) {
                kernel_row_left[a] += weight * left_row[a];
              }
            }
          }
          const float* left_out = left_data +
            static_cast<std::size_t>(row) *
            static_cast<std::size_t>(p_left);
          const float* right_out = symmetric ?
            right_data + static_cast<std::size_t>(row) *
              static_cast<std::size_t>(p_right) :
            NULL;
          for (int b = 0; b < p_right; ++b) {
            const float z = kernel_row_right[b];
            float* output_column =
              output.data() + static_cast<std::size_t>(p_left) * b;
            for (int a = 0; a < p_left; ++a) {
              output_column[a] += left_out[a] * z;
              if (symmetric) {
                output_column[a] += kernel_row_left[a] * right_out[b];
              }
            }
          }
        }
      } catch (...) {
        failed.store(true);
      }
    });
  }
  for (std::thread& worker : workers) worker.join();
  if (failed.load()) stop("A float32 sparse worker failed.");

  NumericMatrix result(p_left, p_right);
  for (int b = 0; b < p_right; ++b) {
    for (int a = 0; a < p_left; ++a) {
      float value = 0.0f;
      const std::size_t position =
        a + static_cast<std::size_t>(p_left) * b;
      for (int thread_index = 0; thread_index < threads; ++thread_index) {
        value += partial[thread_index][position];
      }
      result(a, b) = static_cast<double>(value);
    }
  }
  return result;
}


// Row-parallel K X for downstream kernel consumers. The output is an ordinary
// R double matrix, but products and row accumulation use float32.
//
// [[Rcpp::export(rng = false)]]
NumericMatrix float32_csr_matmul_cpp(
    IntegerVector p,
    IntegerVector j,
    RawVector x,
    IntegerVector dims,
    NumericMatrix input,
    int n_threads = 1,
    bool transpose = false,
    bool symmetric = false) {
  validate_csr(p, j, x, dims, symmetric);
  const int n_rows = dims[0];
  const int n_columns = dims[1];
  const int expected_input_rows = transpose ? n_rows : n_columns;
  if (input.nrow() != expected_input_rows) {
    stop("Input row count does not match the kernel column count.");
  }
  const int n_rhs = input.ncol();

  const int* row_pointer = INTEGER(p);
  const int* column_index = INTEGER(j);
  const Rbyte* values = RAW(x);
  const double* input_pointer = REAL(input);

  if (transpose || symmetric) {
    const int threads = normalized_thread_count(
      n_threads, std::max(1, n_rhs)
    );
    const int result_rows = symmetric ? n_rows : n_columns;
    NumericMatrix result(result_rows, n_rhs);
    double* output_pointer = REAL(result);
    std::vector<std::thread> workers;
    workers.reserve(threads);
    // Each worker allocates its own scratch row, so each can throw bad_alloc.
    // An exception that escapes a std::thread calls std::terminate() and takes
    // the R session with it; catch it and re-raise on the main thread after the
    // join, where stop() is legal.
    std::atomic<bool> worker_failed(false);
    for (int thread_index = 0; thread_index < threads; ++thread_index) {
      const int rhs_begin = static_cast<int>(
        static_cast<std::int64_t>(n_rhs) * thread_index / threads
      );
      const int rhs_end = static_cast<int>(
        static_cast<std::int64_t>(n_rhs) * (thread_index + 1) / threads
      );
      workers.emplace_back([=, &worker_failed]() {
       try {
        std::vector<float> scratch(result_rows, 0.0f);
        for (int rhs = rhs_begin; rhs < rhs_end; ++rhs) {
          std::fill(scratch.begin(), scratch.end(), 0.0f);
          for (int row = 0; row < n_rows; ++row) {
            const float input_value = static_cast<float>(
              input_pointer[
                row + static_cast<std::size_t>(n_rows) * rhs
              ]
            );
            for (int position = row_pointer[row];
                 position < row_pointer[row + 1]; ++position) {
              const int column = column_index[position];
              const float weight = read_float(values, position);
              scratch[column] += weight * input_value;
              if (symmetric) {
                scratch[row] += weight * static_cast<float>(
                  input_pointer[
                    column + static_cast<std::size_t>(n_columns) * rhs
                  ]
                );
              }
            }
          }
          for (int column = 0; column < result_rows; ++column) {
            output_pointer[
              column + static_cast<std::size_t>(result_rows) * rhs
            ] = static_cast<double>(scratch[column]);
          }
        }
       } catch (...) {
        worker_failed.store(true);
       }
      });
    }
    for (std::thread& worker : workers) worker.join();
    if (worker_failed.load()) {
      stop("A float32 sparse matrix-multiply worker failed.");
    }
    return result;
  }

  const int threads = normalized_thread_count(n_threads, std::max(1, n_rows));
  NumericMatrix result(n_rows, n_rhs);
  double* output_pointer = REAL(result);
  std::vector<std::thread> workers;
  workers.reserve(threads);
  std::atomic<bool> worker_failed(false);

  // The transposed/symmetric branch above splits over right-hand sides, so each
  // worker reads one contiguous column of the input and needs no repacking.
  // This branch fans a single nonzero across every right-hand side, which in
  // column-major order means one cache line per side; pack it row-major first.
  const std::vector<float> input_packed =
    pack_row_major_f32(input_pointer, n_columns, n_rhs);
  const float* input_data = input_packed.data();

  for (int thread_index = 0; thread_index < threads; ++thread_index) {
    const int row_begin = static_cast<int>(
      static_cast<std::int64_t>(n_rows) * thread_index / threads
    );
    const int row_end = static_cast<int>(
      static_cast<std::int64_t>(n_rows) * (thread_index + 1) / threads
    );
    workers.emplace_back([=, &worker_failed]() {
     try {
      // Sweep each row's nonzeros once, reading the float32 value and column
      // index a single time and fanning the product across all right-hand
      // sides, instead of re-reading them once per column. The per-right-hand
      // accumulation order is unchanged, so the result is bit-identical.
      std::vector<float> row_output(n_rhs, 0.0f);
      for (int row = row_begin; row < row_end; ++row) {
        std::fill(row_output.begin(), row_output.end(), 0.0f);
        for (int position = row_pointer[row];
             position < row_pointer[row + 1]; ++position) {
          const float weight = read_float(values, position);
          const float* input_row = input_data +
            static_cast<std::size_t>(column_index[position]) *
            static_cast<std::size_t>(n_rhs);
          for (int rhs = 0; rhs < n_rhs; ++rhs) {
            row_output[rhs] += weight * input_row[rhs];
          }
        }
        for (int rhs = 0; rhs < n_rhs; ++rhs) {
          output_pointer[
            row + static_cast<std::size_t>(n_rows) * rhs
          ] = static_cast<double>(row_output[rhs]);
        }
      }
     } catch (...) {
      worker_failed.store(true);
     }
    });
  }
  for (std::thread& worker : workers) worker.join();
  if (worker_failed.load()) {
    stop("A float32 sparse matrix-multiply worker failed.");
  }
  return result;
}


// [[Rcpp::export(rng = false)]]
List float32_csr_sums_cpp(
    IntegerVector p,
    IntegerVector j,
    RawVector x,
    IntegerVector dims,
    bool symmetric = false) {
  validate_csr(p, j, x, dims, symmetric);
  NumericVector row_sums(dims[0]);
  NumericVector column_sums(dims[1]);
  double sum_squares = 0.0;
  const Rbyte* values = RAW(x);
  for (int row = 0; row < dims[0]; ++row) {
    double sum = 0.0;
    for (int position = p[row]; position < p[row + 1]; ++position) {
      const double value = read_float(values, position);
      sum += value;
      column_sums[j[position]] += value;
      sum_squares += value * value;
      if (symmetric) {
        row_sums[j[position]] += value;
        column_sums[row] += value;
        sum_squares += value * value;
      }
    }
    row_sums[row] += sum;
  }
  return List::create(
    _["rowSums"] = row_sums,
    _["colSums"] = column_sums,
    _["sumSquares"] = sum_squares
  );
}


// Diagnostic conversion used by tests and explicit user requests. Large-data
// workflows should keep the kernel encoded and use the operator functions.
//
// [[Rcpp::export(rng = false)]]
NumericVector float32_csr_values_cpp(RawVector x) {
  if (x.size() % 4 != 0) stop("Invalid float32 value-buffer length.");
  NumericVector values(x.size() / 4);
  const Rbyte* source = RAW(x);
  for (R_xlen_t index = 0; index < values.size(); ++index) {
    values[index] = read_float(source, index);
  }
  return values;
}
