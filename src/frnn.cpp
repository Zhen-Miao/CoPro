#include <Rcpp.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <unordered_map>
#include <vector>

using namespace Rcpp;

namespace {

struct CellKey {
  std::int64_t x;
  std::int64_t y;
  std::int64_t z;

  bool operator==(const CellKey& other) const {
    return x == other.x && y == other.y && z == other.z;
  }
};

struct CellKeyHash {
  std::size_t operator()(const CellKey& key) const {
    std::size_t seed = std::hash<std::int64_t>()(key.x);
    seed ^= std::hash<std::int64_t>()(key.y) + 0x9e3779b9U +
            (seed << 6) + (seed >> 2);
    seed ^= std::hash<std::int64_t>()(key.z) + 0x9e3779b9U +
            (seed << 6) + (seed >> 2);
    return seed;
  }
};

std::int64_t checked_grid_size(double span, double radius) {
  const double value = std::floor(span / radius) + 1.0;
  const double max_index = static_cast<double>(
    std::numeric_limits<std::int64_t>::max()
  );
  if (!R_finite(value) || value < 1.0 || value > max_index) {
    stop("Grid extent is too large for the fixed-radius neighbor search.");
  }
  return static_cast<std::int64_t>(value);
}

CellKey point_cell(const NumericMatrix& coords,
                   int row,
                   const std::vector<double>& origin,
                   const std::vector<std::int64_t>& grid_size,
                   double radius) {
  std::int64_t cell[3] = {0, 0, 0};
  const int dimensions = coords.ncol();
  for (int axis = 0; axis < dimensions; ++axis) {
    double raw = std::floor((coords(row, axis) - origin[axis]) / radius);
    if (!R_finite(raw)) {
      stop("Coordinates must contain only finite values.");
    }
    raw = std::max(0.0, raw);
    raw = std::min(raw, static_cast<double>(grid_size[axis] - 1));
    cell[axis] = static_cast<std::int64_t>(raw);
  }
  return CellKey{cell[0], cell[1], cell[2]};
}

bool inside(std::int64_t value, std::int64_t size) {
  return value >= 0 && value < size;
}

}  // namespace


// Exact fixed-radius neighbor enumeration for 2-D and 3-D coordinates.
// The bucket map is used only for direct lookup; traversal order is explicit,
// so identical inputs produce identical output triplets across runs.
// [[Rcpp::export(rng = false)]]
List frnn_grid_cpp(NumericMatrix A, SEXP B_sexp, double radius) {
  const bool within = Rf_isNull(B_sexp);
  NumericMatrix B = within ? A : as<NumericMatrix>(B_sexp);

  const int n_a = A.nrow();
  const int n_b = B.nrow();
  const int dimensions = A.ncol();

  if (dimensions != 2 && dimensions != 3) {
    stop("Fixed-radius neighbor search supports only 2-D or 3-D coordinates.");
  }
  if (B.ncol() != dimensions) {
    stop("A and B must have the same number of coordinate columns.");
  }
  if (!R_finite(radius) || radius <= 0.0) {
    stop("radius r must be a positive finite number");
  }
  if (n_a == 0 || n_b == 0) {
    return List::create(
      _["i"] = IntegerVector(0),
      _["j"] = IntegerVector(0),
      _["d"] = NumericVector(0)
    );
  }

  std::vector<double> origin(dimensions);
  std::vector<double> upper(dimensions);
  std::vector<std::int64_t> grid_size(dimensions);
  for (int axis = 0; axis < dimensions; ++axis) {
    double lo = A(0, axis);
    double hi = A(0, axis);
    if (!R_finite(lo)) stop("Coordinates must contain only finite values.");
    for (int row = 1; row < n_a; ++row) {
      const double value = A(row, axis);
      if (!R_finite(value)) stop("Coordinates must contain only finite values.");
      lo = std::min(lo, value);
      hi = std::max(hi, value);
    }
    for (int row = 0; row < n_b; ++row) {
      const double value = B(row, axis);
      if (!R_finite(value)) stop("Coordinates must contain only finite values.");
      lo = std::min(lo, value);
      hi = std::max(hi, value);
    }
    origin[axis] = lo;
    upper[axis] = hi;
    grid_size[axis] = checked_grid_size(hi - lo, radius);
  }

  std::vector<CellKey> cells_a;
  cells_a.reserve(n_a);
  for (int row = 0; row < n_a; ++row) {
    cells_a.push_back(point_cell(A, row, origin, grid_size, radius));
  }

  std::unordered_map<CellKey, std::vector<int>, CellKeyHash> buckets;
  buckets.reserve(static_cast<std::size_t>(n_b * 1.3) + 1U);
  for (int row = 0; row < n_b; ++row) {
    const CellKey key = point_cell(B, row, origin, grid_size, radius);
    buckets[key].push_back(row);
  }

  std::vector<int> out_i;
  std::vector<int> out_j;
  std::vector<double> out_d;
  const std::size_t reserve_size = std::min<std::size_t>(
    static_cast<std::size_t>(n_a) * 8U, 1000000U
  );
  out_i.reserve(reserve_size);
  out_j.reserve(reserve_size);
  out_d.reserve(reserve_size);

  const double radius_squared = radius * radius;
  for (int a = 0; a < n_a; ++a) {
    if ((a & 16383) == 0) checkUserInterrupt();
    const CellKey base = cells_a[a];

    for (int dz = dimensions == 3 ? -1 : 0;
         dz <= (dimensions == 3 ? 1 : 0); ++dz) {
      const std::int64_t z = base.z + dz;
      if (dimensions == 3 && !inside(z, grid_size[2])) continue;

      for (int dy = -1; dy <= 1; ++dy) {
        const std::int64_t y = base.y + dy;
        if (!inside(y, grid_size[1])) continue;

        for (int dx = -1; dx <= 1; ++dx) {
          const std::int64_t x = base.x + dx;
          if (!inside(x, grid_size[0])) continue;

          const CellKey query{x, y, z};
          const auto bucket = buckets.find(query);
          if (bucket == buckets.end()) continue;

          for (const int b : bucket->second) {
            if (within && a == b) continue;
            double distance_squared = 0.0;
            for (int axis = 0; axis < dimensions; ++axis) {
              const double delta = A(a, axis) - B(b, axis);
              distance_squared += delta * delta;
            }
            if (distance_squared <= radius_squared) {
              out_i.push_back(a + 1);
              out_j.push_back(b + 1);
              out_d.push_back(std::sqrt(distance_squared));
            }
          }
        }
      }
    }
  }

  return List::create(
    _["i"] = wrap(out_i),
    _["j"] = wrap(out_j),
    _["d"] = wrap(out_d)
  );
}
