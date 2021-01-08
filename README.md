# I-PIA

A simple library for fitting an implicit surface by
[progressive-iterative approximation](http://doi.org/10.1016/j.cagd.2020.101817).

Basic geometry is provided by my [libgeom](http://github.com/salvipeter/libgeom/) library.

The test program also uses my [dual contouring](http://github.com/salvipeter/dual-contouring/)
library and the [TCLAP](http://tclap.sourceforge.net/) command line parser.

## Usage

The API consists of a single class called `IPIA`. There are two constructors:
1. `IPIA(degree, bbox, size)`
1. `IPIA(bases, size)`

In the first version, it is assumed that all 3 B-spline bases have the same degree and uniform knots. `bbox` is an array of two points, specifying the minimal and maximal corners of the bounding box. `size` is an array of 3 integers, corresponding to the number of control points in each of the coordinates. The uniform knot vectors are created in a way that the valid range of the B-spline bases covers the whole bounding box.

In the second version, `degree` and `bbox` are substituted by `bases`, which is an array of 3 B-spline bases (`BSBasis` objects from libgeom). This gives the user more flexibility, but care should be taken that the B-spline bases should be defined over the whole fitting area.

The degrees, knot vectors and control network can be queried, and the `()` operator functions as an evaluator.

The fitting operation itself is invoked by

    IPIA surface(...);
    size_t exit_iteration = surface.fit(samples, small_step, iterations, tolerance);
    
... where `samples` contains point/normal pairs (the normals are assumed to be consistently oriented), and `small_step` is a small distance for creating an offset surface (to avoid the trivial solution). There are two exit conditions: `iterations` is the maximum iteration count, and the fitting is also finished when all control points move a distance smaller than `tolerance`. The return value is the actual number of iterations processed.

See also the header file and the example program for more details.
