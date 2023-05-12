#pragma once

#include <cmath>
#include <exception>
#include <geometry.hh>

namespace ImplicitCurvature {

  using namespace Geometry;

  // API (note: when the gradient is 0, these also return 0 instead of nan)

  double mean(const Vector3D &gradient, const Matrix3x3 &hessian);
  double gaussian(const Vector3D &gradient, const Matrix3x3 &hessian);
  std::pair<double, double> principal(const Vector3D &gradient, const Matrix3x3 &hessian);
  std::pair<Vector3D, Vector3D> directions(const Vector3D &gradient, const Matrix3x3 &hessian);


  ////////////////////
  // Implementation //
  ////////////////////

  // Goldman (2005) - equal to (Fuu + Fvv) / (2 |Fn|)
  double mean(const Vector3D &gradient, const Matrix3x3 &hessian) {
    double g2 = gradient.normSqr();
    if (g2 == 0)
      return 0;
    return (g2 * hessian.trace() - gradient * (hessian * gradient)) / (2 * std::pow(g2, 1.5));
  }

  // Goldman (2005) - equal to (Fuu Fvv - Fuv^2) / Fn^2
  double gaussian(const Vector3D &gradient, const Matrix3x3 &hessian) {
    double g2 = gradient.normSqr();
    if (g2 == 0)
      return 0;
    return (gradient * (hessian.adjugate() * gradient)) / (g2 * g2);
  }

  std::pair<double, double> principal(const Vector3D &gradient, const Matrix3x3 &hessian) {
    double km = mean(gradient, hessian);
    double kg = gaussian(gradient, hessian);
    double d = std::sqrt(std::max(km * km - kg, 0.0));
    return { km - d, km + d };
  }

  // Albin et al. (2017)
  std::pair<Vector3D, Vector3D> directions(const Vector3D &gradient, const Matrix3x3 &hessian) {
    if (gradient.normSqr() == 0)
      return { gradient, gradient };

    // Choose an arbitrary (u,v) system in the tangent plane
    Vector3D n = gradient.normalized();
    size_t i = std::abs(n[0]) < std::abs(n[1])
      ? (std::abs(n[0]) < std::abs(n[2]) ? 0 : 2)
      : (std::abs(n[1]) < std::abs(n[2]) ? 1 : 2);
    Vector3D tmp(0, 0, 0); tmp[i] = 1.0;
    Vector3D u = (n ^ tmp).normalize();
    Vector3D v = (n ^ u).normalize();

    // Eq. (15) in the Appendix
    double Fn = gradient.norm();
    double Fuv = u * (hessian * v);
    double Fuu = u * (hessian * u);
    double Fvv = v * (hessian * v);

    auto [kmin, kmax] = principal(gradient, hessian);

    if (std::abs(kmin * Fn - Fuu) >= std::abs(kmin * Fn - Fvv)) {
      double k1 = kmin, k2 = kmax;
      Vector3D t1 = u * Fuv + v * (k1 * Fn - Fuu);
      Vector3D t2 = v * Fuv + u * (k2 * Fn - Fvv);
      return { -t2.normalize(), t1.normalize() };
    }

    double k1 = kmax, k2 = kmin;
    Vector3D t1 = u * Fuv + v * (k1 * Fn - Fuu);
    Vector3D t2 = v * Fuv + u * (k2 * Fn - Fvv);
    return { t1.normalize(), -t2.normalize() };
  }
}
