//! All tests of this crate.

use super::triangulate;
use nalgebra::Vector2;

#[test]
fn triangulate_polygon_with_collinear_vertices() {
    // A quadrilateral with a vertex very close to a corner, on one of its sides.
    let polygon = vec![
        (0, Vector2::new(0.0, 0.0)),
        (1, Vector2::new(2.7e-8, 0.0)),
        (2, Vector2::new(1.0, 0.0)),
        (3, Vector2::new(0.084, 0.916)),
    ];
    assert_eq!(triangulate(&polygon).map(|t| t.len()), Some(2));
    // Many collinear vertices on one side.
    let polygon = vec![
        (0, Vector2::new(0.0, 0.0)),
        (1, Vector2::new(0.25, 0.0)),
        (2, Vector2::new(0.5, 0.0)),
        (3, Vector2::new(0.75, 0.0)),
        (4, Vector2::new(1.0, 0.0)),
        (5, Vector2::new(0.0, 1.0)),
    ];
    assert_eq!(triangulate(&polygon).map(|t| t.len()), Some(4));
}

#[test]
fn triangulate_triangle_with_collinear_vertices_on_its_sides() {
    let polygon = vec![
        (0, Vector2::new(0.0, 0.0)),
        (1, Vector2::new(0.5, 0.0)),
        (2, Vector2::new(1.0, 0.0)),
        (3, Vector2::new(0.5, 0.5)),
        (4, Vector2::new(0.0, 1.0)),
        (5, Vector2::new(0.0, 0.5)),
    ];
    let triangles = triangulate(&polygon).unwrap();
    assert_eq!(triangles.len(), polygon.len() - 2);
    for &(id, _) in &polygon {
        assert!(triangles.iter().flatten().any(|&v| v == id));
    }
    for [a, b, c] in triangles {
        let (a, b, c) = (polygon[a].1, polygon[b].1, polygon[c].1);
        assert!((b - a).perp(&(c - a)) > 0.);
    }
}

type Polygon = Vec<(usize, Vector2<f64>)>;

fn area2(a: Vector2<f64>, b: Vector2<f64>, c: Vector2<f64>) -> f64 {
    (b - a).perp(&(c - a))
}

/// Ids are offset from the indices, so that the ids (not the indices) must end up in the triangles.
fn polygon(points: &[Vector2<f64>]) -> Polygon {
    points.iter().enumerate().map(|(i, &p)| (1000 + 3 * i, p)).collect()
}

/// Triangulates and checks: n - 2 counter-clockwise triangles of positive area, using exactly the polygon's ids
/// (every one of them), whose total area is the polygon's area.
fn check(polygon: &Polygon) {
    let n = polygon.len();
    let triangles = triangulate(polygon).unwrap_or_else(|| panic!("no triangulation of {polygon:?}"));
    assert_eq!(triangles.len(), n - 2);
    let position = |id: usize| polygon.iter().find(|&&(i, _)| i == id).expect("unknown id").1;
    for &(id, _) in polygon {
        assert!(triangles.iter().flatten().any(|&v| v == id), "id {id} unused");
    }
    let mut total = 0.;
    for &[a, b, c] in &triangles {
        let area = area2(position(a), position(b), position(c));
        assert!(area > 0., "triangle {a} {b} {c} has area {area}");
        total += area;
    }
    let expected: f64 = (0..n).map(|i| polygon[i].1.perp(&polygon[(i + 1) % n].1)).sum();
    assert!(
        (total - expected).abs() <= 1e-9 * expected.abs(),
        "{total} vs {expected}"
    );
}

/// A (counter-clockwise) triangle with the given parameters (in (0, 1), sorted) of points on each side.
fn cut_triangle(corners: [Vector2<f64>; 3], sides: [&[f64]; 3]) -> Polygon {
    let mut points = vec![];
    for (i, side) in sides.iter().enumerate() {
        let (a, b) = (corners[i], corners[(i + 1) % 3]);
        points.push(a);
        points.extend(side.iter().map(|&t| a + (b - a) * t));
    }
    polygon(&points)
}

/// A small deterministic xorshift generator of values in [0, 1).
struct Rng(u64);
impl Rng {
    fn next(&mut self) -> f64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        (self.0 >> 11) as f64 / (1u64 << 53) as f64
    }
    fn sorted(&mut self, count: usize) -> Vec<f64> {
        let mut values: Vec<f64> = (0..count).map(|_| self.next()).filter(|&t| t > 0.).collect();
        values.sort_by(f64::total_cmp);
        values.dedup();
        values
    }
}

const CORNERS: [Vector2<f64>; 3] = [Vector2::new(0., 0.), Vector2::new(1., 0.), Vector2::new(0.084, 0.916)];

#[test]
fn many_collinear_vertices_on_all_sides() {
    let uniform: Vec<f64> = (1..=50).map(|i| f64::from(i) / 51.).collect();
    check(&cut_triangle(CORNERS, [&uniform, &uniform, &uniform]));
    // Very uneven spacing: points 1e-9 apart, and points very close to the corners (as in the mesh cutting).
    let clustered: Vec<f64> = (0..20).map(|i| 0.5 + 1e-9 * f64::from(i)).collect();
    let near_corners = [2.7e-8, 5e-8, 0.3, 1. - 1e-7, 1. - 2.7e-8];
    check(&cut_triangle(CORNERS, [&clustered, &near_corners, &[2.7e-8]]));
    check(&cut_triangle(CORNERS, [&near_corners, &near_corners, &clustered]));
}

/// Rounding makes the points on a side slightly non-collinear (turns of about +-1e-26), so the polygon is not exactly
/// convex. This used to fail (tiny ears were clipped until only collinear vertices remained); the exact turn tests now
/// keep the remainder from degenerating.
#[test]
fn rounded_collinear_vertices_close_together() {
    let side = [0.500_000_015, 0.500_000_016, 0.500_000_017, 0.500_000_018];
    check(&cut_triangle(CORNERS, [&[3e-9], &[0.999_999_999], &side]));
}

#[test]
fn sliver_triangle_with_vertices_on_its_sides() {
    let sliver = [Vector2::new(0., 0.), Vector2::new(1., 0.), Vector2::new(0.5, 1e-6)];
    check(&cut_triangle(sliver, [&[0.1, 0.5, 0.9], &[0.5], &[0.25, 0.75]]));
    check(&cut_triangle(sliver, [&[0.5], &[], &[]]));
    check(&cut_triangle(sliver, [&[], &[2.7e-8, 0.5], &[0.5, 1. - 2.7e-8]]));
}

#[test]
fn regular_and_large_convex_polygons() {
    let regular = |n: usize| {
        let points: Vec<_> = (0..n)
            .map(|i| {
                let angle = std::f64::consts::TAU * i as f64 / n as f64;
                Vector2::new(angle.cos(), angle.sin())
            })
            .collect();
        polygon(&points)
    };
    for n in 3..=100 {
        check(&regular(n));
    }
    check(&regular(1000));
}

#[test]
fn random_convex_polygons_and_cut_triangles() {
    let mut rng = Rng(0x9e37_79b9_7f4a_7c15);
    for _ in 0..2000 {
        // Points on an ellipse at random angles.
        let (rx, ry) = (0.01 + rng.next(), 0.01 + rng.next());
        let count = 3 + (rng.next() * 30.) as usize;
        let points: Vec<_> = rng
            .sorted(count)
            .iter()
            .map(|t| {
                Vector2::new(
                    rx * (std::f64::consts::TAU * t).cos(),
                    ry * (std::f64::consts::TAU * t).sin(),
                )
            })
            .collect();
        if points.len() >= 3 {
            check(&polygon(&points));
        }
    }
    for _ in 0..2000 {
        let corners = [0, 1, 2].map(|_| Vector2::new(rng.next(), rng.next()));
        if area2(corners[0], corners[1], corners[2]) < 1e-6 {
            continue;
        }
        let mut side = || {
            let count = (rng.next() * 8.) as usize;
            rng.sorted(count)
        };
        let sides = [side(), side(), side()];
        check(&cut_triangle(corners, [&sides[0], &sides[1], &sides[2]]));
    }
}

#[test]
fn vertices_on_one_side_and_minimal_cases() {
    let fan: Vec<f64> = (1..100).map(|i| f64::from(i) / 100.).collect();
    for sides in [[&fan[..], &[], &[]], [&[], &fan[..], &[]], [&[], &[], &fan[..]]] {
        check(&cut_triangle(CORNERS, sides));
    }
    check(&cut_triangle(CORNERS, [&[], &[], &[]]));
    check(&cut_triangle(CORNERS, [&[0.5], &[], &[]]));
    check(&polygon(&[
        Vector2::new(0., 0.),
        Vector2::new(1., 0.),
        Vector2::new(1., 1.),
        Vector2::new(0., 1.),
    ]));
}

#[test]
fn degenerate_polygons_have_no_triangulation() {
    let collinear = |n: usize| {
        polygon(
            &(0..n)
                .map(|i| Vector2::new(i as f64, 2. * i as f64))
                .collect::<Vec<_>>(),
        )
    };
    for n in 3..10 {
        assert_eq!(triangulate(&collinear(n)), None);
    }
    // Clockwise.
    let clockwise = polygon(&[Vector2::new(0., 0.), Vector2::new(0., 1.), Vector2::new(1., 0.)]);
    assert_eq!(triangulate(&clockwise), None);
}

#[test]
fn non_convex_polygons() {
    let v = Vector2::new;
    // An L, with collinear points on its long sides.
    check(&polygon(&[
        v(0., 0.),
        v(1., 0.),
        v(2., 0.),
        v(2., 1.),
        v(1., 1.),
        v(1., 2.),
        v(1., 3.),
        v(0., 3.),
        v(0., 1.5),
    ]));
    // A comb: teeth between deep notches.
    let mut comb = vec![v(0., 0.), v(7., 0.), v(7., 3.)];
    for k in (0..7).rev() {
        let x = f64::from(k);
        comb.extend([v(x + 0.5, 1.), v(x, 3.)]);
    }
    check(&polygon(&comb));
}

#[test]
fn random_star_shaped_polygons() {
    // Vertices at random radii around a center, in angular order: simple, and mostly not convex.
    let mut rng = Rng(7);
    for n in [5, 8, 13, 30, 100] {
        for _ in 0..50 {
            let points = (0..n)
                .map(|i| {
                    let angle = std::f64::consts::TAU * f64::from(i) / f64::from(n);
                    let radius = 0.2 + rng.next();
                    Vector2::new(radius * angle.cos(), radius * angle.sin())
                })
                .collect::<Vec<_>>();
            check(&polygon(&points));
        }
    }
}

#[test]
fn thin_triangle_whose_area_cancels_from_the_far_corner() {
    let v = Vector2::new;
    let far = v(0.610_942_326_830_717_7, 0.920_930_634_801_503_6);
    let (a, b) = (
        v(0.215_264_216_664_856_6, 0.246_064_354_324_545_1),
        v(0.215_264_216_664_856_7, 0.246_064_354_324_545_24),
    );
    assert!(triangulate(&polygon(&[far, a, b])).is_some());
}
