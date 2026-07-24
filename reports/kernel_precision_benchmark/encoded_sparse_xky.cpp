#include <Rcpp.h>
#include <cstdint>
#include <cstring>
#include <vector>

using namespace Rcpp;

namespace {

uint32_t as_bits(float value) {
  uint32_t bits;
  std::memcpy(&bits, &value, sizeof(bits));
  return bits;
}

float from_bits(uint32_t bits) {
  float value;
  std::memcpy(&value, &bits, sizeof(value));
  return value;
}

// IEEE-754 binary32 -> binary16, round-to-nearest-even.
uint16_t float_to_half(float value) {
  const uint32_t bits = as_bits(value);
  const uint32_t sign = (bits >> 16) & 0x8000u;
  const uint32_t exponent = (bits >> 23) & 0xffu;
  uint32_t mantissa = bits & 0x7fffffu;

  if (exponent == 0xffu) {
    if (mantissa == 0u) return static_cast<uint16_t>(sign | 0x7c00u);
    mantissa >>= 13;
    return static_cast<uint16_t>(
      sign | 0x7c00u | mantissa | (mantissa == 0u)
    );
  }

  const int half_exponent = static_cast<int>(exponent) - 127 + 15;
  if (half_exponent >= 31) {
    return static_cast<uint16_t>(sign | 0x7c00u);
  }

  if (half_exponent <= 0) {
    if (half_exponent < -10) return static_cast<uint16_t>(sign);
    mantissa |= 0x800000u;
    const int shift = 14 - half_exponent;
    uint32_t half_mantissa = mantissa >> shift;
    const uint32_t remainder = mantissa & ((1u << shift) - 1u);
    const uint32_t halfway = 1u << (shift - 1);
    if (remainder > halfway ||
        (remainder == halfway && (half_mantissa & 1u))) {
      ++half_mantissa;
    }
    return static_cast<uint16_t>(sign | half_mantissa);
  }

  uint32_t half = sign |
    (static_cast<uint32_t>(half_exponent) << 10) |
    (mantissa >> 13);
  const uint32_t remainder = mantissa & 0x1fffu;
  if (remainder > 0x1000u || (remainder == 0x1000u && (half & 1u))) {
    ++half;
  }
  return static_cast<uint16_t>(half);
}

float half_to_float(uint16_t half) {
  const uint32_t sign = (static_cast<uint32_t>(half & 0x8000u)) << 16;
  uint32_t exponent = (half >> 10) & 0x1fu;
  uint32_t mantissa = half & 0x3ffu;
  uint32_t bits;

  if (exponent == 0u) {
    if (mantissa == 0u) {
      bits = sign;
    } else {
      int e = -14;
      while ((mantissa & 0x400u) == 0u) {
        mantissa <<= 1;
        --e;
      }
      mantissa &= 0x3ffu;
      bits = sign |
        (static_cast<uint32_t>(e + 127) << 23) |
        (mantissa << 13);
    }
  } else if (exponent == 0x1fu) {
    bits = sign | 0x7f800000u | (mantissa << 13);
  } else {
    bits = sign | ((exponent - 15u + 127u) << 23) | (mantissa << 13);
  }
  return from_bits(bits);
}

uint16_t read_u16(const RawVector& bytes, R_xlen_t index) {
  const R_xlen_t offset = index * 2;
  return static_cast<uint16_t>(bytes[offset]) |
    (static_cast<uint16_t>(bytes[offset + 1]) << 8);
}

uint32_t read_u32(const RawVector& bytes, R_xlen_t index) {
  const R_xlen_t offset = index * 4;
  return static_cast<uint32_t>(bytes[offset]) |
    (static_cast<uint32_t>(bytes[offset + 1]) << 8) |
    (static_cast<uint32_t>(bytes[offset + 2]) << 16) |
    (static_cast<uint32_t>(bytes[offset + 3]) << 24);
}

float decode_value(
    const RawVector& values,
    R_xlen_t index,
    int encoding,
    const NumericVector& codebook) {
  if (encoding == 1) {
    return half_to_float(read_u16(values, index));
  }
  if (encoding == 2) {
    return from_bits(read_u32(values, index));
  }
  const unsigned int code = values[index];
  return static_cast<float>(codebook[code]);
}

void validate_sparse_inputs(
    const IntegerVector& p,
    const IntegerVector& i,
    const RawVector& values,
    const IntegerVector& dims,
    const NumericMatrix& x_left,
    const NumericMatrix& x_right,
    int encoding,
    const NumericVector& codebook) {
  if (dims.size() != 2) stop("dims must have length two");
  if (p.size() != dims[1] + 1) stop("invalid compressed-column pointer");
  if (x_left.nrow() != dims[0]) stop("x_left row count does not match K");
  if (x_right.nrow() != dims[1]) stop("x_right row count does not match K");
  if (i.size() != p[p.size() - 1]) stop("invalid sparse row-index length");
  const R_xlen_t expected = encoding == 1 ? i.size() * 2 :
    (encoding == 2 ? i.size() * 4 : i.size());
  if (values.size() != expected) stop("encoded value length does not match nnz");
  if (encoding == 3 && codebook.size() == 0) stop("bin codebook is empty");
}

NumericMatrix encoded_xky_impl(
    const IntegerVector& p,
    const IntegerVector& i,
    const RawVector& values,
    const IntegerVector& dims,
    const NumericMatrix& x_left,
    const NumericMatrix& x_right,
    int encoding,
    const NumericVector& codebook,
    bool float_accumulation) {
  validate_sparse_inputs(
    p, i, values, dims, x_left, x_right, encoding, codebook
  );

  const int n_rows = dims[0];
  const int n_cols = dims[1];
  const int n_left = x_left.ncol();
  const int n_right = x_right.ncol();
  NumericMatrix result(n_left, n_right);

  if (float_accumulation) {
    std::vector<float> kx(
      static_cast<size_t>(n_rows) * static_cast<size_t>(n_right), 0.0f
    );
    for (int col = 0; col < n_cols; ++col) {
      for (int pos = p[col]; pos < p[col + 1]; ++pos) {
        const int row = i[pos];
        const float weight = decode_value(values, pos, encoding, codebook);
        for (int b = 0; b < n_right; ++b) {
          kx[row + static_cast<size_t>(n_rows) * b] +=
            weight * static_cast<float>(x_right(col, b));
        }
      }
    }

    std::vector<float> y(
      static_cast<size_t>(n_left) * static_cast<size_t>(n_right), 0.0f
    );
    for (int b = 0; b < n_right; ++b) {
      for (int row = 0; row < n_rows; ++row) {
        const float z = kx[row + static_cast<size_t>(n_rows) * b];
        for (int a = 0; a < n_left; ++a) {
          y[a + static_cast<size_t>(n_left) * b] +=
            static_cast<float>(x_left(row, a)) * z;
        }
      }
    }
    for (int b = 0; b < n_right; ++b) {
      for (int a = 0; a < n_left; ++a) {
        result(a, b) = y[a + static_cast<size_t>(n_left) * b];
      }
    }
    return result;
  }

  std::vector<double> kx(
    static_cast<size_t>(n_rows) * static_cast<size_t>(n_right), 0.0
  );
  for (int col = 0; col < n_cols; ++col) {
    for (int pos = p[col]; pos < p[col + 1]; ++pos) {
      const int row = i[pos];
      const double weight = decode_value(values, pos, encoding, codebook);
      for (int b = 0; b < n_right; ++b) {
        kx[row + static_cast<size_t>(n_rows) * b] +=
          weight * x_right(col, b);
      }
    }
  }

  for (int b = 0; b < n_right; ++b) {
    for (int row = 0; row < n_rows; ++row) {
      const double z = kx[row + static_cast<size_t>(n_rows) * b];
      for (int a = 0; a < n_left; ++a) {
        result(a, b) += x_left(row, a) * z;
      }
    }
  }
  return result;
}

}  // namespace

// [[Rcpp::export]]
RawVector cpp_pack_half(NumericVector values) {
  RawVector bytes(values.size() * 2);
  for (R_xlen_t index = 0; index < values.size(); ++index) {
    const uint16_t half = float_to_half(static_cast<float>(values[index]));
    bytes[index * 2] = half & 0xffu;
    bytes[index * 2 + 1] = (half >> 8) & 0xffu;
  }
  return bytes;
}

// [[Rcpp::export]]
NumericVector cpp_unpack_half(RawVector bytes) {
  if (bytes.size() % 2 != 0) stop("half byte vector has odd length");
  NumericVector values(bytes.size() / 2);
  for (R_xlen_t index = 0; index < values.size(); ++index) {
    values[index] = half_to_float(read_u16(bytes, index));
  }
  return values;
}

// [[Rcpp::export]]
RawVector cpp_pack_float32(NumericVector values) {
  RawVector bytes(values.size() * 4);
  for (R_xlen_t index = 0; index < values.size(); ++index) {
    const uint32_t bits = as_bits(static_cast<float>(values[index]));
    bytes[index * 4] = bits & 0xffu;
    bytes[index * 4 + 1] = (bits >> 8) & 0xffu;
    bytes[index * 4 + 2] = (bits >> 16) & 0xffu;
    bytes[index * 4 + 3] = (bits >> 24) & 0xffu;
  }
  return bytes;
}

// [[Rcpp::export]]
NumericVector cpp_unpack_float32(RawVector bytes) {
  if (bytes.size() % 4 != 0) stop("float32 byte vector has invalid length");
  NumericVector values(bytes.size() / 4);
  for (R_xlen_t index = 0; index < values.size(); ++index) {
    values[index] = from_bits(read_u32(bytes, index));
  }
  return values;
}

// [[Rcpp::export]]
NumericMatrix cpp_xky_half(
    IntegerVector p,
    IntegerVector i,
    RawVector values,
    IntegerVector dims,
    NumericMatrix x_left,
    NumericMatrix x_right,
    bool float_accumulation = false) {
  return encoded_xky_impl(
    p, i, values, dims, x_left, x_right, 1, NumericVector(),
    float_accumulation
  );
}

// [[Rcpp::export]]
NumericMatrix cpp_xky_float32(
    IntegerVector p,
    IntegerVector i,
    RawVector values,
    IntegerVector dims,
    NumericMatrix x_left,
    NumericMatrix x_right,
    bool float_accumulation = false) {
  return encoded_xky_impl(
    p, i, values, dims, x_left, x_right, 2, NumericVector(),
    float_accumulation
  );
}

// [[Rcpp::export]]
NumericMatrix cpp_xky_bins(
    IntegerVector p,
    IntegerVector i,
    RawVector codes,
    NumericVector codebook,
    IntegerVector dims,
    NumericMatrix x_left,
    NumericMatrix x_right,
    bool float_accumulation = false) {
  return encoded_xky_impl(
    p, i, codes, dims, x_left, x_right, 3, codebook,
    float_accumulation
  );
}
