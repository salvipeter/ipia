#include "ipia.hh"

#include <Eigen/Sparse>

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

double infinityNorm(const Eigen::SparseMatrix<double> &m) {
  double result = 0;
  for (size_t i = 0; i < (size_t)m.rows(); ++i) {
    double sum = 0;
    for (size_t j = 0; j < (size_t)m.cols(); ++j)
      sum += std::abs(m.coeff(i, j));
    if (sum > result)
      result = sum;
  }
  return result;
}

void IPIA::fit(const std::vector<PointNormal> &samples, double small_step, size_t iterations) {
  // Variable names follow the notations of Section 3.2
  // Note that here we set epsilon = sigma.
  cpts.assign(sizes[0] * sizes[1] * sizes[2], 0.0);
  double sigma = small_step, e = sigma;
  size_t n = samples.size(), N = cpts.size();
  // Set up the collocation matrix
  std::vector<Eigen::Triplet<double>> triplets;
  auto addTriplets = [&](size_t row, const Point3D &p) {
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
          triplets.emplace_back(row, index, coeff[0][i] * coeff[1][j] * coeff[2][k]);
        }
      }
    }
  };
  for (size_t pi = 0; pi < n; ++pi) {
    addTriplets(pi, samples[pi].p);
    addTriplets(n + pi, samples[pi].p + samples[pi].n * sigma);
  }
  Eigen::SparseMatrix<double> B(2 * n, N);
  B.setFromTriplets(triplets.begin(), triplets.end());
  Eigen::SparseMatrix<double> Bt = B.transpose();
  Eigen::SparseMatrix<double> BtB = Bt * B;
  Eigen::MatrixXd I = Eigen::MatrixXd::Identity(N, N);
  Eigen::Map<Eigen::VectorXd> C(&cpts[0], N);
  double mu = 2.0 / infinityNorm(BtB);
  Eigen::VectorXd b = Eigen::VectorXd::Zero(2 * n);
  b.tail(n) = Eigen::VectorXd::Constant(n, e);
  for (size_t i = 0; i < iterations; ++i)
    C = (I - mu * BtB) * C + mu * Bt * b;
}
