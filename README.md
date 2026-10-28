# vangogh

Triangulation of simple polygons by ear clipping, robust to collinear vertices. Given a simple polygon, convex or
not, whose boundary may contain arbitrarily many collinear vertices (such as a mesh triangle subdivided along its
edges), it computes a triangulation that uses every vertex and contains no degenerate (zero-area) triangles.

<img src="gallery.png" alt="Triangulations of polygons with collinear vertices: a cut triangle, a fan, a hexagon, an L, a comb and a star" width="600">

## Usage

```rust
use nalgebra::Vector2;

// A counter-clockwise polygon of (id, position): a triangle with two extra vertices on its bottom side.
let polygon = [
    (0, Vector2::new(0., 0.)),
    (1, Vector2::new(1., 0.)),
    (2, Vector2::new(2., 0.)),
    (3, Vector2::new(3., 0.)),
    (4, Vector2::new(0., 3.)),
];
let triangles = vangogh::triangulate(&polygon).unwrap(); // three triangles of ids, counter-clockwise
```

`triangulate` returns `None` only if the input is not a simple counter-clockwise polygon with positive area, e.g. when
all vertices are collinear or the order is clockwise. A planar 3D polygon can be triangulated in a projection onto its
plane; the ids carry the triangles back.

## How it works

This crate builds on two results:

- **Two ears** guarantee that ear clipping succeeds: Gary H. Meisters. *[Polygons Have
  Ears](https://doi.org/10.2307/2319703)*. The American Mathematical Monthly 82(6), 1975.
- **Exact predicates** decide every orientation without rounding errors: Jonathan Richard Shewchuk. *[Adaptive
  Precision Floating-Point Arithmetic and Fast Robust Geometric Predicates](https://doi.org/10.1007/PL00009321)*.
  Discrete & Computational Geometry 18(3), 1997.

An ear is clipped only if no other vertex lies in it and both the ear and the rest of the polygon keep a positive area.
Among the valid ears, the one that leaves the largest smaller area goes first, which avoids slivers.

```bibtex
@article{Meisters:1975:PHE,
  author  = {Meisters, Gary H.},
  title   = {Polygons Have Ears},
  journal = {The American Mathematical Monthly}, volume = {82}, number = {6}, pages = {648--651}, year = {1975}
}
@article{Shewchuk:1997:APF,
  author  = {Shewchuk, Jonathan Richard},
  title   = {Adaptive Precision Floating-Point Arithmetic and Fast Robust Geometric Predicates},
  journal = {Discrete \& Computational Geometry}, volume = {18}, number = {3}, pages = {305--363}, year = {1997}
}
```

## License

MIT
