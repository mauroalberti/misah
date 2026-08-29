
use super::plane::Plane;
use super::point::Point3D;
use super::segment::Segment3D;

/// A triangle in space, as its three vertices.
///
/// Vertex order fixes the winding, and so the sign of the normal that
/// `to_plane` returns; nothing here depends on it, and neither does the
/// attitude read back by `GeologicalPlane::from_plane`, which orients the
/// normal upward regardless.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Triangle3D {
    pub vertices: [Point3D; 3],
}

impl Triangle3D {

    pub fn new(a: Point3D, b: Point3D, c: Point3D) -> Self {
        Self { vertices: [a, b, c] }
    }

    pub fn vertex(&self, i: usize) -> Option<Point3D> {
        self.vertices.get(i).copied()
    }

    pub fn area(&self) -> f64 {

        let v1 = self.vertices[0].vector_to(&self.vertices[1]);
        let v2 = self.vertices[0].vector_to(&self.vertices[2]);

        v1.cross(&v2).norm() / 2.0
    }

    /// The plane through the three vertices.
    ///
    /// `None` when the triangle is degenerate -- two vertices coincident, or
    /// the three collinear. `Plane::from_three_points` decides that on the sine
    /// of the angle between the two edges, which is scale-invariant: a fixed
    /// cutoff on the area instead (as the C++ original used, at 1e-10) means
    /// something quite different for a triangle a metre across than for one
    /// spanning a UTM sheet.
    pub fn to_plane(&self) -> Option<Plane> {
        Plane::from_three_points(self.vertices[0], self.vertices[1], self.vertices[2])
    }

    /// The three sides, each starting where the previous one ends.
    pub fn edges(&self) -> [Segment3D; 3] {
        [
            Segment3D::new(self.vertices[0], self.vertices[1]),
            Segment3D::new(self.vertices[1], self.vertices[2]),
            Segment3D::new(self.vertices[2], self.vertices[0]),
        ]
    }

    /// The axis-aligned bounding box, as its lower and upper corners.
    pub fn bounds(&self) -> (Point3D, Point3D) {

        let mut min = self.vertices[0].coords;
        let mut max = min;

        for v in &self.vertices[1..] {
            for k in 0..3 {
                if v.coords[k] < min[k] { min[k] = v.coords[k]; }
                if v.coords[k] > max[k] { max[k] = v.coords[k]; }
            }
        }

        (Point3D::from(min), Point3D::from(max))
    }

    /// Whether a point already known to lie in the triangle's plane falls
    /// inside the triangle, boundary included.
    ///
    /// The test is the sign of three cross products, so it says nothing about
    /// how far off the plane the point may be: a point far above the triangle
    /// but within its prism answers `true`. Callers get coplanarity from how
    /// the point was built -- in the mesh-grid kernel, by interpolating to the
    /// plane along a DEM edge -- rather than from a tolerance here, which at
    /// UTM coordinates would have to be chosen against the wrong scale.
    pub fn contains_coplanar_point(&self, point: &Point3D) -> bool {

        let [a, b, c] = self.vertices;

        same_side(point, &a, &b, &c) && same_side(point, &b, &a, &c) && same_side(point, &c, &a, &b)
    }
}

/// Whether `p1` and `p2` lie on the same side of the line through `a` and `b`.
///
/// Both cross products are taken against the same edge vector, so comparing
/// their directions answers the question without normalizing either.
fn same_side(p1: &Point3D, p2: &Point3D, a: &Point3D, b: &Point3D) -> bool {

    let edge = a.vector_to(b);

    let cross_1 = edge.cross(&a.vector_to(p1));
    let cross_2 = edge.cross(&a.vector_to(p2));

    cross_1.dot(&cross_2) >= 0.0
}

#[cfg(test)]
mod tests {
    use super::*;

    fn unit_triangle() -> Triangle3D {
        Triangle3D::new(
            Point3D::from([0.0, 0.0, 0.0]),
            Point3D::from([1.0, 0.0, 0.0]),
            Point3D::from([0.0, 1.0, 0.0]),
        )
    }

    #[test]
    fn area_of_the_unit_right_triangle_is_one_half() {
        assert!((unit_triangle().area() - 0.5).abs() < 1e-12);
    }

    #[test]
    fn area_of_a_collapsed_triangle_is_zero() {
        let t = Triangle3D::new(
            Point3D::from([0.0, 0.0, 0.0]),
            Point3D::from([1.0, 0.0, 0.0]),
            Point3D::from([2.0, 0.0, 0.0]),
        );

        assert_eq!(t.area(), 0.0);
    }

    #[test]
    fn to_plane_passes_through_every_vertex() {
        let t = Triangle3D::new(
            Point3D::from([0.0, 0.0, 1.0]),
            Point3D::from([3.0, 0.0, 2.0]),
            Point3D::from([0.0, 4.0, -1.0]),
        );
        let plane = t.to_plane().expect("a well-conditioned triangle");

        for v in &t.vertices {
            assert!(plane.distance_to_point(v) < 1e-12);
        }
    }

    #[test]
    fn to_plane_of_a_collinear_triangle_is_none() {
        let t = Triangle3D::new(
            Point3D::from([0.0, 0.0, 0.0]),
            Point3D::from([1.0, 1.0, 1.0]),
            Point3D::from([2.0, 2.0, 2.0]),
        );

        assert!(t.to_plane().is_none());
    }

    #[test]
    fn edges_form_a_closed_ring() {
        let t = unit_triangle();
        let e = t.edges();

        assert_eq!(e[0].end_pt, e[1].start_pt);
        assert_eq!(e[1].end_pt, e[2].start_pt);
        assert_eq!(e[2].end_pt, e[0].start_pt);
    }

    #[test]
    fn bounds_span_all_three_vertices() {
        let t = Triangle3D::new(
            Point3D::from([1.0, -2.0, 5.0]),
            Point3D::from([-3.0, 4.0, 0.0]),
            Point3D::from([2.0, 1.0, 7.0]),
        );

        let (min, max) = t.bounds();

        assert_eq!(min, Point3D::from([-3.0, -2.0, 0.0]));
        assert_eq!(max, Point3D::from([2.0, 4.0, 7.0]));
    }

    #[test]
    fn an_interior_point_is_contained() {
        assert!(unit_triangle().contains_coplanar_point(&Point3D::from([0.25, 0.25, 0.0])));
    }

    #[test]
    fn a_point_outside_the_triangle_is_not_contained() {
        assert!(!unit_triangle().contains_coplanar_point(&Point3D::from([1.0, 1.0, 0.0])));
    }

    #[test]
    fn a_point_on_an_edge_counts_as_contained() {
        // The boundary is inclusive: an intersection point landing exactly on a
        // mesh triangle's edge belongs to it, rather than falling through the
        // gap between two adjacent triangles.
        assert!(unit_triangle().contains_coplanar_point(&Point3D::from([0.5, 0.0, 0.0])));
    }

    #[test]
    fn a_vertex_counts_as_contained() {
        assert!(unit_triangle().contains_coplanar_point(&Point3D::from([0.0, 1.0, 0.0])));
    }

    #[test]
    fn containment_holds_for_a_tilted_triangle() {
        // Same test, off the coordinate planes, so nothing can pass by having
        // dropped a coordinate.
        let t = Triangle3D::new(
            Point3D::from([0.0, 0.0, 0.0]),
            Point3D::from([2.0, 0.0, 2.0]),
            Point3D::from([0.0, 2.0, 2.0]),
        );

        assert!(t.contains_coplanar_point(&Point3D::from([0.5, 0.5, 1.0])));
        assert!(!t.contains_coplanar_point(&Point3D::from([2.0, 2.0, 4.0])));
    }
}
