//! Ear clipping triangulation of a simple polygon whose vertices may be collinear, such as a mesh triangle with extra
//! vertices on its sides where it is cut by curves.
//!
//! Input: a simple polygon (convex or not) in counter-clockwise order; consecutive collinear vertices are allowed. Not
//! for polygons with holes or self-intersections. A 3D polygon in a plane can be triangulated in a projection onto that
//! plane (e.g., dropping the coordinate of its normal's largest component); the ids carry the triangles back.
//!
//! Output: `n - 2` triangles (counter-clockwise) that use every vertex and have positive area, or `None` if no such
//! triangulation is found (e.g., all vertices are collinear). An ear is only clipped if no other vertex lies in it and
//! both the ear and the remaining polygon have positive area; among those, our heuristic picks the ear that maximizes
//! the smaller of the two areas, which avoids slivers.
//!
//! Whether an ear or the remaining polygon has positive area is decided exactly, with Shewchuk's adaptive precision
//! orientation predicate ([1997](https://doi.org/10.1007/PL00009321)), so rounding in nearly collinear vertices cannot
//! leave a flat remainder.
//!
//! Ear clipping relies on the two-ears theorem: G. H. Meisters, "Polygons have ears", American Mathematical Monthly
//! 82(6), 648-651, 1975, <https://doi.org/10.2307/2319703>.
//!
//! ```
//! use nalgebra::Vector2;
//!
//! // A triangle with two extra vertices on its bottom side.
//! let polygon = [
//!     (0, Vector2::new(0., 0.)),
//!     (1, Vector2::new(1., 0.)),
//!     (2, Vector2::new(2., 0.)),
//!     (3, Vector2::new(3., 0.)),
//!     (4, Vector2::new(0., 3.)),
//! ];
//! let triangles = vangogh::triangulate(&polygon).unwrap();
//! assert_eq!(triangles.len(), 3);
//! ```

use nalgebra::Vector2;
use ordered_float::OrderedFloat;
use std::collections::BinaryHeap;

type Vector2D = Vector2<f64>;

// Positive if `a`, `b` and `c` turn counter-clockwise, zero if they are collinear (exactly).
fn orient(a: Vector2D, b: Vector2D, c: Vector2D) -> f64 {
    let coord = |v: Vector2D| robust::Coord { x: v.x, y: v.y };
    robust::orient2d(coord(a), coord(b), coord(c))
}

// The corner opposite the longest side: computing the area from it is most accurate in floating point.
fn base(t: [Vector2D; 3]) -> usize {
    let sides = [0, 1, 2].map(|k| (t[(k + 2) % 3] - t[(k + 1) % 3]).norm_squared());
    (0..3).max_by(|&i, &j| sides[i].total_cmp(&sides[j])).unwrap()
}

// Twice the area of a triangle, from its base corner.
fn area(mut t: [Vector2D; 3]) -> f64 {
    let k = base(t);
    t.rotate_left(k);
    (t[1] - t[0]).perp(&(t[2] - t[0]))
}

/// Triangulate a counter-clockwise polygon of (id, position). None if that fails (e.g., all vertices are collinear).
///
/// Greedy: with exact predicates, a simple polygon always has a valid ear, and what remains after clipping it is again
/// a simple polygon with positive area (two-ears theorem), so no choice leads to a dead end.
pub fn triangulate(polygon: &[(usize, Vector2D)]) -> Option<Vec<[usize; 3]>> {
    let mut polygon = polygon.to_vec();
    let indices = |n: usize, i: usize| [(i + n - 1) % n, i, (i + 1) % n];
    let corners = |polygon: &[(usize, Vector2D)], i: usize| indices(polygon.len(), i).map(|k| polygon[k].1);
    // The triangle at `i`, starting at its base corner (so its area is positive in floating point from its first corner).
    let triangle = |polygon: &[(usize, Vector2D)], i: usize| {
        let mut ids = indices(polygon.len(), i).map(|k| polygon[k].0);
        ids.rotate_left(base(corners(polygon, i)));
        ids
    };
    let turn = |polygon: &[(usize, Vector2D)], i: usize| {
        let [a, b, c] = corners(polygon, i);
        orient(a, b, c)
    };
    let n = polygon.len();
    // Per corner: its exact turn and its area (cached: clipping an ear only changes those of its neighbors).
    let mut turns: Vec<f64> = (0..n).map(|i| turn(&polygon, i)).collect();
    let mut areas: Vec<f64> = (0..n).map(|i| area(corners(&polygon, i))).collect();
    let mut total = (0..n).map(|i| polygon[i].1.perp(&polygon[(i + 1) % n].1)).sum::<f64>() / 2.;
    let mut triangles = Vec::with_capacity(n.saturating_sub(2));
    while polygon.len() > 3 {
        let n = polygon.len();
        // Whether no other vertex lies in the ear (or on its sides). Only non-convex vertices need testing: if any vertex
        // lies in the ear, one of those does too. Outside the bounding box is (exactly) outside.
        let empty = |i: usize| {
            let [a, b, c] = corners(&polygon, i);
            let (min, max) = (a.inf(&b).inf(&c), a.sup(&b).sup(&c));
            (0..n)
                .filter(|&k| turns[k] <= 0. && (k + n + 1 - i) % n > 2)
                .map(|k| polygon[k].1)
                .all(|q| {
                    q != q.sup(&min).inf(&max) || orient(a, b, q) < 0. || orient(b, c, q) < 0. || orient(c, a, q) < 0.
                })
        };
        // The best ear: the smaller of its and the rest's area is largest (ties: last index), and it has a normal.
        let mut candidates: BinaryHeap<_> = (0..n)
            .filter(|&i| turns[i] > 0. && areas[i] > 0.)
            .map(|i| (OrderedFloat(areas[i].min(total - areas[i] / 2.)), i))
            .collect();
        let ear = std::iter::from_fn(|| candidates.pop())
            .map(|(_, i)| i)
            .find(|&i| empty(i))?;
        total -= areas[ear] / 2.;
        triangles.push(triangle(&polygon, ear));
        polygon.remove(ear);
        turns.remove(ear);
        areas.remove(ear);
        for k in [(ear + n - 2) % (n - 1), ear % (n - 1)] {
            (turns[k], areas[k]) = (turn(&polygon, k), area(corners(&polygon, k)));
        }
    }
    let ok = polygon.len() == 3 && turns[1] > 0. && areas[1] > 0.;
    ok.then(|| triangles.push(triangle(&polygon, 1)));
    ok.then_some(triangles)
}

#[cfg(test)]
mod tests;
