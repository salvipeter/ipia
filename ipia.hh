#pragma once

// Based on:
// Y.F. Hamza, H. Lin, Z. Li:
//   Implicit progressive-iterative approximation for curve and surface reconstruction.
//     Computer Aided Geometric Design 77, #101817, 2020.
// https://doi.org/10.1016/j.cagd.2020.101817

#include <functional>

#include <geometry.hh>

class IPIA {
public:
  // Constructors
  IPIA(size_t degree, const std::array<Geometry::Point3D, 2> &bbox,
       const std::array<size_t, 3> &size); // equal degrees, uniform knots
  IPIA(const std::array<Geometry::BSBasis, 3> &bases, const std::array<size_t, 3> &size);
  IPIA(const IPIA &) = default;
  IPIA &operator=(const IPIA &) = default;

  // Properties
  size_t size(size_t i) const;                         // # of control points  in the i-th axis
  size_t degree(size_t i) const;                       //   B-spline degree    in the i-th axis
  const Geometry::DoubleVector &knots(size_t i) const; // B-spline knot vector in the i-th axis
  double controlPoint(size_t i, size_t j, size_t k) const;
  const Geometry::DoubleVector &controlNet() const;    // row-major

  // Evaluation
  double operator()(const Geometry::Point3D &p) const;
  double derivative(const Geometry::Point3D &p, size_t dx, size_t dy, size_t dz) const;
  Geometry::Vector3D gradient(const Geometry::Point3D &p) const;
  Geometry::Matrix3x3 hessian(const Geometry::Point3D &p) const;

  // Fitting
  struct PointNormal {
    Geometry::Point3D p;
    Geometry::Vector3D n;
  };
  size_t fit(const std::vector<PointNormal> &samples, double small_step,
             size_t iterations, double tolerance);

private:
  void doBasis(const Geometry::Point3D &p, const std::function<void(size_t, double)> &f) const;
  void doBasisDerivative(const Geometry::Point3D &p, const std::function<void(size_t, double)> &f,
                         size_t dx, size_t dy, size_t dz) const;
  double computeMu(const Geometry::PointVector &points);

  std::array<size_t, 3> sizes;
  std::array<Geometry::BSBasis, 3> bases;
  Geometry::DoubleVector cpts;
};
