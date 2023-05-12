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

size_t IPIA::size(size_t i) const {
  return sizes[i];
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

void IPIA::doBasis(const Point3D &p, const std::function<void(size_t, double)> &f) const {
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
        double B = coeff[0][i] * coeff[1][j] * coeff[2][k];
        f(index, B);
      }
    }
  }
}

void IPIA::doBasisDerivative(const Point3D &p, const std::function<void(size_t, double)> &f,
                             size_t dx, size_t dy, size_t dz) const {
  size_t degree[3], span[3], d[3] = { dx, dy, dz };
  DoubleMatrix coeff[3];
  for (size_t i = 0; i < 3; ++i) {
    degree[i] = bases[i].degree();
    span[i] = bases[i].findSpan(p[i]);
    bases[i].basisFunctionDerivatives(span[i], p[i], d[i], coeff[i]);
  }
  for (size_t i = 0; i <= degree[0]; ++i) {
    size_t x_offset = (span[0] - degree[0] + i) * sizes[1] * sizes[2];
    for (size_t j = 0; j <= degree[1]; ++j) {
      size_t y_offset = (span[1] - degree[1] + j) * sizes[2];
      for (size_t k = 0; k <= degree[2]; ++k) {
        size_t z_offset = span[2] - degree[2] + k;
        size_t index = x_offset + y_offset + z_offset;
        double B = coeff[0][d[0]][i] * coeff[1][d[1]][j] * coeff[2][d[2]][k];
        f(index, B);
      }
    }
  }
}

double IPIA::operator()(const Point3D &p) const {
  double result = 0;
  doBasis(p, [&](size_t index, double B) { result += cpts[index] * B; });
  return result;
}

double IPIA::derivative(const Geometry::Point3D &p, size_t dx, size_t dy, size_t dz) const {
  double result = 0;
  doBasisDerivative(p, [&](size_t index, double B) { result += cpts[index] * B; }, dx, dy, dz);
  return result;
}

Vector3D IPIA::gradient(const Geometry::Point3D &p) const {
  return { derivative(p, 1, 0, 0), derivative(p, 0, 1, 0), derivative(p, 0, 0, 1) };
}

Matrix3x3 IPIA::hessian(const Geometry::Point3D &p) const {
  double xx = derivative(p, 2, 0, 0);
  double xy = derivative(p, 1, 1, 0);
  double xz = derivative(p, 1, 0, 1);
  double yy = derivative(p, 0, 2, 0);
  double yz = derivative(p, 0, 1, 1);
  double zz = derivative(p, 0, 0, 2);
  return Matrix3x3({ xx, xy, xz, xy, yy, yz, xz, yz, zz });
}

double IPIA::computeMu(const PointVector &points) {
  DoubleVector sums(sizes[0] * sizes[1] * sizes[2], 0);
  for (const auto &p : points)
    doBasis(p, [&](size_t index, double B) { sums[index] += B; });
  return 2.0 / *std::max_element(sums.begin(), sums.end());
}

size_t IPIA::fit(const std::vector<PointNormal> &samples, double small_step,
                 size_t iterations, double tolerance) {
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

  size_t iter = 0;
  for (; iter < iterations; ++iter) {
    double max_delta = 0;
    for (size_t pi = 0; pi < 2 * n; ++pi)
      delta[pi] = (pi < n ? 0 : e) - (*this)(points[pi]);
    for (size_t pi = 0; pi < 2 * n; ++pi)
      doBasis(points[pi], [&](size_t index, double B) { Delta[index] += B * delta[pi]; });
    for (size_t i = 0; i < N; ++i) {
      if (std::abs(Delta[i]) > max_delta)
        max_delta = std::abs(Delta[i]);
      cpts[i] += mu * Delta[i];
      Delta[i] = 0;
    }
    if (max_delta < tolerance)
      break;
  }
  return iter;
}
