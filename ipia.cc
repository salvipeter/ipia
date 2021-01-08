#include "ipia.hh"

#include <algorithm>

using namespace Geometry;

IPIA::IPIA(size_t degree, const std::array<Point3D, 2> &bbox, const std::array<size_t, 3> &size) {
  for (size_t i = 0; i < 3; ++i) {
    sizes[i] = size[i];
    bases[i].setDegree(degree);
    auto &knots = bases[i].knots();
    // Compute the knot range such that evaluation is possible in the whole bounding box
    double len = bbox[1][i] - bbox[0][i];
    double extra = len * degree / (sizes[i] - degree);
    double min = bbox[0][i] - extra;
    len += extra * 2;
    for (size_t j = 0; j <= sizes[i] + degree; ++j) {
      double u = (double)j / (sizes[i] + degree);
      knots.push_back(min + len * u);
    }
  }
}

IPIA::IPIA(const std::array<BSBasis, 3> &bases, const std::array<size_t, 3> &size)
  : sizes(size), bases(bases) {
}


size_t IPIA::degree(size_t i) const {
  return bases[i].degree();
}

const DoubleVector &IPIA::knots(size_t i) const {
  return bases[i].knots();
}

double IPIA::controlPoint(size_t i, size_t j, size_t k) const {
  return cpts[i*sizes[1]*sizes[2]+j*sizes[2]+k];
}

const DoubleVector &IPIA::controlNet() const {
  return cpts;
}

double IPIA::operator()(const Point3D &p) const {
  size_t degree[3], span[3];
  DoubleVector coeff[3];
  double result = 0;
  for (size_t i = 0; i < 3; ++i) {
    degree[i] = bases[i].degree();
    span[i] = bases[i].findSpan(p[i]);
    bases[i].basisFunctions(span[i], p[i], coeff[i]);
  }
  for (size_t i = 0; i <= degree[0]; ++i) {
    size_t x_offset = (span[0] - degree[0] + i) * sizes[1] * sizes[2];
    for (size_t j = 0; j <= degree[1]; ++j) {
      size_t y_offset = (span[1] - degree[1] + j) * sizes[2];
      for (size_t k = 0; k <= degree[2]; ++k) {
        size_t z_offset = span[2] - degree[2] + k;
        size_t index = x_offset + y_offset + z_offset;
        result += cpts[index] * coeff[0][i] * coeff[1][j] * coeff[2][k];
      }
    }
  }
  return result;
}

double IPIA::computeMu(const PointVector &points) {
  DoubleVector sums(sizes[0] * sizes[1] * sizes[2], 0);
  for (const auto &p : points) {
    size_t degree[3], span[3];
    DoubleVector coeff[3];
    for (size_t i = 0; i < 3; ++i) {
      degree[i] = bases[i].degree();
      span[i] = bases[i].findSpan(p[i]);
      bases[i].basisFunctions(span[i], p[i], coeff[i]);
    }
    for (size_t i = 0; i <= degree[0]; ++i) {
      size_t x_offset = (span[0] - degree[0] + i) * sizes[1] * sizes[2];
      for (size_t j = 0; j <= degree[1]; ++j) {
        size_t y_offset = (span[1] - degree[1] + j) * sizes[2];
        for (size_t k = 0; k <= degree[2]; ++k) {
          size_t z_offset = span[2] - degree[2] + k;
          size_t index = x_offset + y_offset + z_offset;
          sums[index] += coeff[0][i] * coeff[1][j] * coeff[2][k];
        }
      }
    }
  }
  return 2.0 / *std::max_element(sums.begin(), sums.end());
}

void IPIA::fit(const std::vector<PointNormal> &samples, double small_step, size_t iterations) {
  // Variable names follow the notations of Section 3.2
  // Note that here we set epsilon = sigma.
  cpts.assign(sizes[0] * sizes[1] * sizes[2], 0.0);
  const size_t n = samples.size(), N = cpts.size();
  const double sigma = small_step, e = sigma;
  PointVector points(2 * n);
  for (size_t i = 0; i < n; ++i) {
    points[i] = samples[i].p;
    points[n+i] = samples[i].p + samples[i].n * sigma;
  };
  const double mu = computeMu(points);
  DoubleVector delta(2 * n), Delta(N, 0.0);

  for (size_t iter = 0; iter < iterations; ++iter) {
    for (size_t pi = 0; pi < 2 * n; ++pi)
      delta[pi] = (pi < n ? 0 : e) - (*this)(points[pi]);
    for (size_t pi = 0; pi < 2 * n; ++pi) {
      const auto &p = points[pi];
      size_t degree[3], span[3];
      DoubleVector coeff[3];
      for (size_t i = 0; i < 3; ++i) {
        degree[i] = bases[i].degree();
        span[i] = bases[i].findSpan(p[i]);
        bases[i].basisFunctions(span[i], p[i], coeff[i]);
      }
      for (size_t i = 0; i <= degree[0]; ++i) {
        size_t x_offset = (span[0] - degree[0] + i) * sizes[1] * sizes[2];
        for (size_t j = 0; j <= degree[1]; ++j) {
          size_t y_offset = (span[1] - degree[1] + j) * sizes[2];
          for (size_t k = 0; k <= degree[2]; ++k) {
            size_t z_offset = span[2] - degree[2] + k;
            size_t index = x_offset + y_offset + z_offset;
            Delta[index] += coeff[0][i] * coeff[1][j] * coeff[2][k] * delta[pi];
          }
        }
      }
    }
    for (size_t i = 0; i < N; ++i) {
      cpts[i] += mu * Delta[i];
      Delta[i] = 0;
    }
  }
}
