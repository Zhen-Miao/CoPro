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

void write_float(Rbyte* destination, R_xlen_t index, float value) {
  std::uint32_t bits;
  std::memcpy(&bits, &value, sizeof(bits));
  const R_xlen_t offset = index * 4;
  destination[offset] = bits & 0xffu;
  destination[offset + 1] = (bits >> 8) & 0xffu;
  destination[offset + 2] = (bits >> 16) & 0xffu;
  destination[offset + 3] = (bits >> 24) & 0xffu;
}

float read_float(const Rbyte* source, R_xlen_t index) {
  const R_xlen_t offset = index * 4;
  const std::uint32_t bits =
    static_cast<std::uint32_t>(source[offset]) |
    (static_cast<std::uint32_t>(source[offset + 1]) << 8) |
    (static_cast<std::uint32_t>(source[offset + 2]) << 16) |
    (static_cast<std::uint32_t>(source[offset + 3]) << 24);
  float value;
  std::memcpy(&value, &bits, sizeof(value));
  return value;
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

template <typename Callback>
std::size_t visit_neighbors(
    int row_begin,
    int row_end,
    const NumericMatrix& A,
    const NumericMatrix& B,
    const std::vector<CellKeyF32>& cells_a,
    const std::unordered_map<
      CellKeyF32, std::vector<int>, CellKeyF32Hash
    >& buckets,
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
          const auto bucket = buckets.find(CellKeyF32{x, y, z});
          if (bucket == buckets.end()) continue;
          for (const int b : bucket->second) {
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

void validate_csr(
    const IntegerVector& p,
    const IntegerVector& j,
    const RawVector& x,
    const IntegerVector& dims,
    bool symmetric = false) {
  if (dims.size() != 2) stop("dims must have length two.");
  if (dims[0] < 0 || dims[1] < 0) stop("dims must be nonnegative.");
  if (p.size() != dims[0] + 1) {
    stop("CSR row pointer length does not match nrow.");
  }
  if (p[0] != 0 || p[p.size() - 1] != j.size()) {
    stop("Invalid CSR row pointer.");
  }
  if (x.size() != static_cast<R_xlen_t>(j.size()) * 4) {
    stop("Float32 value buffer does not match the CSR nonzero count.");
  }
  if (symmetric && dims[0] != dims[1]) {
    stop("A symmetric CSR kernel must be square.");
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
    int normalization = 0) {
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
  std::unordered_map<
    CellKeyF32, std::vector<int>, CellKeyF32Hash
  > buckets;
  buckets.reserve(static_cast<std::size_t>(n_b * 1.3) + 1U);
  for (int row = 0; row < n_b; ++row) {
    buckets[
      point_cell_f32(B, row, origin, grid_size, radius)
    ].push_back(row);
  }

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
    stop("A single float32 sparse kernel block would exceed the 32-bit ",
         "compressed-index limit.");
  }

  std::vector<FloatEdge> edges;
  edges.reserve(static_cast<std::size_t>(reserve_double));
  std::vector<int> edge_pointer(static_cast<std::size_t>(n_a) + 1U, 0);
  std::size_t zero_distance_count = 0U;
  float minimum_nonzero = std::numeric_limits<float>::infinity();

  for (int row = 0; row < n_a; ++row) {
    if ((row & 16383) == 0) checkUserInterrupt();
    visit_neighbors(
      row, row + 1, A, B, cells_a, buckets, grid_size,
      radius_squared, percentile, scaling_factor,
      truncate_low_distance, symmetric,
      [&](int, int column, float distance) {
        if (distance == 0.0f) {
          ++zero_distance_count;
        } else {
          minimum_nonzero = std::min(minimum_nonzero, distance);
        }
        edges.push_back(FloatEdge{column, distance});
      }
    );
    if (edges.size() >
        static_cast<std::size_t>(std::numeric_limits<int>::max())) {
      stop("A single float32 sparse kernel block exceeds the 32-bit ",
           "compressed-index limit.");
    }
    edge_pointer[static_cast<std::size_t>(row) + 1U] =
      static_cast<int>(edges.size());
  }
  if (zero_distance_count > 0U &&
      std::isfinite(static_cast<double>(minimum_nonzero))) {
    for (FloatEdge& edge : edges) {
      if (edge.distance == 0.0f) edge.distance = minimum_nonzero;
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
      for (int position = edge_pointer[row];
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
        for (int position = edge_pointer[row];
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
        for (int position = edge_pointer[row];
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
          stop("A normalized float32 sparse kernel exceeds the 32-bit ",
               "index limit.");
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
        for (int position = edge_pointer[row];
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
        for (int position = edge_pointer[row];
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
    _["temporary_bytes"] =
      static_cast<double>(edges.capacity() * sizeof(FloatEdge)) +
      static_cast<double>(edge_pointer.capacity() * sizeof(int)),
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
  const double* left = REAL(x_left);
  const double* right = REAL(x_right);

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
            for (int b = 0; b < p_right; ++b) {
              kernel_row_right[b] += weight * static_cast<float>(
                right[column + static_cast<std::size_t>(n_columns) * b]
              );
            }
            if (symmetric) {
              for (int a = 0; a < p_left; ++a) {
                kernel_row_left[a] += weight * static_cast<float>(
                  left[column + static_cast<std::size_t>(n_rows) * a]
                );
              }
            }
          }
          for (int b = 0; b < p_right; ++b) {
            const float z = kernel_row_right[b];
            for (int a = 0; a < p_left; ++a) {
              output[
                a + static_cast<std::size_t>(p_left) * b
              ] += static_cast<float>(
                left[row + static_cast<std::size_t>(n_rows) * a]
              ) * z;
              if (symmetric) {
                output[
                  a + static_cast<std::size_t>(p_left) * b
                ] += kernel_row_left[a] * static_cast<float>(
                  right[
                    row + static_cast<std::size_t>(n_columns) * b
                  ]
                );
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
    for (int thread_index = 0; thread_index < threads; ++thread_index) {
      const int rhs_begin = static_cast<int>(
        static_cast<std::int64_t>(n_rhs) * thread_index / threads
      );
      const int rhs_end = static_cast<int>(
        static_cast<std::int64_t>(n_rhs) * (thread_index + 1) / threads
      );
      workers.emplace_back([=]() {
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
      });
    }
    for (std::thread& worker : workers) worker.join();
    return result;
  }

  const int threads = normalized_thread_count(n_threads, std::max(1, n_rows));
  NumericMatrix result(n_rows, n_rhs);
  double* output_pointer = REAL(result);
  std::vector<std::thread> workers;
  workers.reserve(threads);

  for (int thread_index = 0; thread_index < threads; ++thread_index) {
    const int row_begin = static_cast<int>(
      static_cast<std::int64_t>(n_rows) * thread_index / threads
    );
    const int row_end = static_cast<int>(
      static_cast<std::int64_t>(n_rows) * (thread_index + 1) / threads
    );
    workers.emplace_back([=]() {
      for (int row = row_begin; row < row_end; ++row) {
        for (int rhs = 0; rhs < n_rhs; ++rhs) {
          float value = 0.0f;
          for (int position = row_pointer[row];
               position < row_pointer[row + 1]; ++position) {
            value += read_float(values, position) * static_cast<float>(
              input_pointer[
                column_index[position] +
                  static_cast<std::size_t>(n_columns) * rhs
              ]
            );
          }
          output_pointer[
            row + static_cast<std::size_t>(n_rows) * rhs
          ] = static_cast<double>(value);
        }
      }
    });
  }
  for (std::thread& worker : workers) worker.join();
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
