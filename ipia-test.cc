#include <chrono>
#include <cmath>
#include <fstream>

#include <tclap/CmdLine.h>

#include "dc.hh"
#include "implicit-curvatures.hh"
#include "ipia.hh"

using namespace Geometry;

std::array<Point3D, 2> boundingBox(const TriMesh &mesh) {
  Point3D boxmin, boxmax;
  const auto &points = mesh.points();
  boxmin = boxmax = points[0];
  for (const auto &p : points)
    for (int i = 0; i < 3; ++i) {
      boxmin[i] = std::min(boxmin[i], p[i]);
      boxmax[i] = std::max(boxmax[i], p[i]);
    }
  // Add 5%
  auto mean = (boxmin + boxmax) / 2;
  boxmin = mean + (boxmin - mean) * 1.05;
  boxmax = mean + (boxmax - mean) * 1.05;
  return { boxmin, boxmax };
}

std::array<size_t, 3> computeResolution(const std::array<Point3D, 2> &bbox, size_t size) {
  std::array<size_t, 3> resolution;
  auto axis = bbox[1] - bbox[0];
  double axis_delta = axis.norm() / size / std::sqrt(3);
  resolution[0] = std::max<size_t>((size_t)std::ceil(axis[0] / axis_delta), 2);
  resolution[1] = std::max<size_t>((size_t)std::ceil(axis[1] / axis_delta), 2);
  resolution[2] = std::max<size_t>((size_t)std::ceil(axis[2] / axis_delta), 2);
  return resolution;
}

void approximateNormals(std::vector<IPIA::PointNormal> &pns,
                        const std::list<TriMesh::Triangle> &tris) {
  // Weights according to:
  //   N. Max, Weights for computing vertex normals from facet normals.
  //     Journal of Graphics Tools, Vol. 4(2), 1999.
  for (const auto &tri : tris) {
    for (size_t i = 0; i < 3; ++i) {
      size_t i0 = tri[i], i1 = tri[(i+1)%3], i2 = tri[(i+2)%3];
      Vector3D v1 = pns[i0].p - pns[i2].p, v2 = pns[i1].p - pns[i0].p;
      double w = v1.normSqr() * v2.normSqr();
      pns[i0].n += (v1 ^ v2) / (w == 0.0 ? 1.0 : w);
    }
  }
  for (auto &pn : pns)
    if (pn.n.norm() > epsilon)
      pn.n.normalize();
}

void writeVertexCurvatures(const std::vector<Point3D> &vertices, const IPIA &surface,
                           std::string filename) {
  std::ofstream f(filename);
  f << "# vtk DataFile Version 2.0" << std::endl;
  f << "Vertices with normals, mean, Gaussian and principal curvature values & directions"
    << std::endl;
  f << "ASCII" << std::endl;
  f << "DATASET POLYDATA" << std::endl;
  f << "POINTS " << vertices.size() << " float" << std::endl;
  for (const auto &v : vertices)
    f << v << std::endl;
  f << "POINT_DATA " << vertices.size() << std::endl;
  f << "NORMALS normal float" << std::endl;
  for (const auto &v : vertices) {
    auto g = surface.gradient(v);
    if (g.normSqr() > 0)
      g.normalize();
    f << g << std::endl;
  }
  f << "SCALARS mean float 1" << std::endl;
  f << "LOOKUP_TABLE default" << std::endl;
  for (const auto &v : vertices)
    f << ImplicitCurvature::mean(surface.gradient(v), surface.hessian(v)) << std::endl;
  f << "SCALARS Gaussian float 1" << std::endl;
  f << "LOOKUP_TABLE default" << std::endl;
  for (const auto &v : vertices)
    f << ImplicitCurvature::gaussian(surface.gradient(v), surface.hessian(v)) << std::endl;
  f << "SCALARS k1 float 1" << std::endl;
  f << "LOOKUP_TABLE default" << std::endl;
  for (const auto &v : vertices)
    f << ImplicitCurvature::principal(surface.gradient(v), surface.hessian(v)).first << std::endl;
  f << "SCALARS k2 float 1" << std::endl;
  f << "LOOKUP_TABLE default" << std::endl;
  for (const auto &v : vertices)
    f << ImplicitCurvature::principal(surface.gradient(v), surface.hessian(v)).second << std::endl;
  f << "NORMALS d1 float" << std::endl;
  for (const auto &v : vertices)
    f << ImplicitCurvature::directions(surface.gradient(v), surface.hessian(v)).first << std::endl;
  f << "NORMALS d2 float" << std::endl;
  for (const auto &v : vertices)
    f << ImplicitCurvature::directions(surface.gradient(v), surface.hessian(v)).second << std::endl;
}

int main(int argc, char **argv) {
  TCLAP::CmdLine cmd(R"END(Implicit PIA fitting test
Given an input mesh, it approximates normal vectors, and creates an implicit surface via PIA fitting. Then the output quadmesh is created by dual contouring.)END", ' ', "0.1");
  TCLAP::ValueArg<std::string> infileArg("i", "input", "Input mesh", true, "", "input.obj", cmd);
  TCLAP::ValueArg<std::string> outfileArg("o", "output", "Output filename", false, "/tmp/output",
                                          "output", cmd);
  TCLAP::ValueArg<size_t> controlArg("c", "control", "B-spline control size", false, 10,
                                     "# of CPs", cmd);
  TCLAP::ValueArg<size_t> resArg("r", "resolution", "Output resolution", false, 50,
                                 "# of cells", cmd);
  TCLAP::ValueArg<size_t> iterArg("n", "iterations", "Max. iteration count", false, 10,
                                  "# of iterations", cmd);
  TCLAP::ValueArg<double> stepArg("s", "step", "Step size", false, 0.01,
                                  "distance", cmd);
  TCLAP::ValueArg<double> tolArg("t", "tolerance", "Maximum deviation at convergence", false, 1e-5,
                                  "deviation", cmd);

  try {
    cmd.parse(argc, argv);
  } catch (TCLAP::ArgException &e) {
    std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
    return 1;
  }

  size_t control = controlArg.getValue(), res = resArg.getValue(), iterations = iterArg.getValue();
  double step = stepArg.getValue(), tol = tolArg.getValue();
  std::string infile = infileArg.getValue(), outfile = outfileArg.getValue();

  auto mesh = TriMesh::readOBJ(infile);
  std::cout << mesh.points().size() << " points loaded." << std::endl;

  std::chrono::steady_clock::time_point start, stop;

  auto bbox = boundingBox(mesh);
  auto size = computeResolution(bbox, control);
  auto dc_res = computeResolution(bbox, res);
  std::cout << "Bounding box: (" << bbox[0] << ") - (" << bbox[1] << ')' << std::endl;
  std::cout << "Control size: " << size[0] << 'x' << size[1] << 'x' << size[2]
            << " (= " << size[0] * size[1] * size[2] << " CPs)" << std::endl;
  std::cout << "Dual contouring resolution: "
            << dc_res[0] << 'x' << dc_res[1] << 'x' << dc_res[2]
            << " (= " << dc_res[0] * dc_res[1] * dc_res[2] << " cells)" << std::endl;

  start = std::chrono::steady_clock::now();
  std::vector<IPIA::PointNormal> samples;
  size_t n = mesh.points().size();
  samples.resize(n);
  for (size_t i = 0; i < n; ++i)
    samples[i] = { mesh[i], { 0, 0, 0 } };
  approximateNormals(samples, mesh.triangles());
  stop = std::chrono::steady_clock::now();
  std::cout << "Creating input: "
            << std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count()
            << "ms" << std::endl;

  start = std::chrono::steady_clock::now();
  IPIA surface(3, bbox, size);
  size_t done_iterations = surface.fit(samples, step, iterations, tol);
  stop = std::chrono::steady_clock::now();
  std::cout << "Fitting: "
            << std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count()
            << "ms" << std::endl;
  if (done_iterations == iterations)
    std::cout << "Reached maximum iteration count" << std::endl;
  else
    std::cout << "Converged at " << done_iterations << " iterations" << std::endl;

  start = std::chrono::steady_clock::now();
  auto eval = [&](const DualContouring::Point3D &p) { return surface({ p[0], p[1], p[2] }); };
  std::array<DualContouring::Point3D, 2> dc_bbox = { {
    { bbox[0][0], bbox[0][1], bbox[0][2] },
    { bbox[1][0], bbox[1][1], bbox[1][2] }
  } };
  auto dc_mesh = DualContouring::isosurface(eval, 0.0, dc_bbox, dc_res);
  dc_mesh.writeOBJ(outfile + ".obj");
  stop = std::chrono::steady_clock::now();
  std::cout << "Mesh generation: "
            << std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count()
            << "ms" << std::endl;

  std::cout << "Mesh written to " << outfile + ".obj" << std::endl;

  start = std::chrono::steady_clock::now();
  PointVector vertices;
  std::transform(dc_mesh.points.begin(), dc_mesh.points.end(), std::back_inserter(vertices),
                 [](const auto &p) { return Point3D(p[0], p[1], p[2]); });
  writeVertexCurvatures(vertices, surface, outfile + ".vtk");  
  stop = std::chrono::steady_clock::now();
  std::cout << "Curvature generation: "
            << std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count()
            << "ms" << std::endl;

  std::cout << "Curvatures written to " << outfile + ".vtk" << std::endl;
}
