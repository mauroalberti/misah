
use thiserror::Error;

use super::point::Point3D;
use super::triangle::Triangle3D;

#[derive(Debug, Error)]
pub enum MeshError {
    #[error("vertex index {index} is out of range for a mesh of {vertices} vertices")]
    VertexOutOfRange { index: usize, vertices: usize },
}

/// A triangular surface, as a shared vertex pool plus the index triplets that
/// form the faces.
///
/// Indexed rather than a flat list of `Triangle3D` because that is how meshes
/// arrive -- a VTK `POLYDATA` file gives points once and refers to them by
/// index -- and because the vertex sharing is the reason a mesh face count is
/// roughly twice its vertex count rather than three times.
#[derive(Debug, Clone, Default)]
pub struct TriangleMesh {
    vertices: Vec<Point3D>,
    faces: Vec<[usize; 3]>,
}

impl TriangleMesh {

    /// Build a mesh from vertices and face indices, checking every index.
    ///
    /// The check is done once here so the accessors need not repeat it, and so
    /// a malformed file fails at the point it is read rather than somewhere
    /// downstream in a kernel.
    pub fn new(vertices: Vec<Point3D>, faces: Vec<[usize; 3]>) -> Result<Self, MeshError> {

        for face in &faces {
            for &index in face {
                if index >= vertices.len() {
                    return Err(MeshError::VertexOutOfRange {
                        index,
                        vertices: vertices.len(),
                    });
                }
            }
        }

        Ok(Self { vertices, faces })
    }

    /// Build a mesh from VTK-style triangle strips.
    ///
    /// Within a strip, every three consecutive indices form a triangle, so a
    /// strip of n points yields n-2 of them; strips shorter than three points
    /// contribute none. The alternating winding a strip implies is not applied:
    /// consecutive triangles therefore face opposite ways, which is immaterial
    /// here since attitudes are read off an upward-oriented normal, and it
    /// keeps the faces in step with the C++ original this was ported from.
    pub fn from_triangle_strips(
        vertices: Vec<Point3D>,
        strips: &[Vec<usize>],
    ) -> Result<Self, MeshError> {

        let mut faces = Vec::new();

        for strip in strips {
            for window in strip.windows(3) {
                faces.push([window[0], window[1], window[2]]);
            }
        }

        Self::new(vertices, faces)
    }

    /// A mesh of one face, for a caller holding a single plane patch rather
    /// than a surface.
    pub fn from_triangle(triangle: &Triangle3D) -> Self {
        Self {
            vertices: triangle.vertices.to_vec(),
            faces: vec![[0, 1, 2]],
        }
    }

    pub fn vertices(&self) -> &[Point3D] {
        &self.vertices
    }

    pub fn faces(&self) -> &[[usize; 3]] {
        &self.faces
    }

    pub fn num_triangles(&self) -> usize {
        self.faces.len()
    }

    pub fn is_empty(&self) -> bool {
        self.faces.is_empty()
    }

    pub fn triangle(&self, index: usize) -> Option<Triangle3D> {

        let [a, b, c] = *self.faces.get(index)?;

        Some(Triangle3D::new(self.vertices[a], self.vertices[b], self.vertices[c]))
    }

    pub fn triangles(&self) -> impl Iterator<Item = Triangle3D> + '_ {
        self.faces
            .iter()
            .map(move |&[a, b, c]| {
                Triangle3D::new(self.vertices[a], self.vertices[b], self.vertices[c])
            })
    }

    /// The axis-aligned bounding box of the vertex pool, or `None` when empty.
    ///
    /// Taken over all vertices, including any not referenced by a face.
    pub fn bounds(&self) -> Option<(Point3D, Point3D)> {

        let first = self.vertices.first()?;

        let mut min = first.coords;
        let mut max = min;

        for v in &self.vertices[1..] {
            for k in 0..3 {
                if v.coords[k] < min[k] { min[k] = v.coords[k]; }
                if v.coords[k] > max[k] { max[k] = v.coords[k]; }
            }
        }

        Some((Point3D::from(min), Point3D::from(max)))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn square_vertices() -> Vec<Point3D> {
        vec![
            Point3D::from([0.0, 0.0, 0.0]),
            Point3D::from([1.0, 0.0, 0.0]),
            Point3D::from([0.0, 1.0, 0.0]),
            Point3D::from([1.0, 1.0, 0.0]),
        ]
    }

    #[test]
    fn a_face_index_past_the_vertex_pool_is_rejected() {
        let err = TriangleMesh::new(square_vertices(), vec![[0, 1, 9]]).unwrap_err();

        assert!(matches!(
            err,
            MeshError::VertexOutOfRange { index: 9, vertices: 4 }
        ));
    }

    #[test]
    fn a_strip_of_n_points_yields_n_minus_two_triangles() {
        let mesh =
            TriangleMesh::from_triangle_strips(square_vertices(), &[vec![0, 1, 2, 3]]).unwrap();

        assert_eq!(mesh.num_triangles(), 2);
        assert_eq!(mesh.faces(), &[[0, 1, 2], [1, 2, 3]]);
    }

    #[test]
    fn strips_shorter_than_a_triangle_contribute_nothing() {
        let mesh =
            TriangleMesh::from_triangle_strips(square_vertices(), &[vec![0, 1], vec![2], vec![]])
                .unwrap();

        assert!(mesh.is_empty());
    }

    #[test]
    fn several_strips_are_concatenated() {
        let mesh = TriangleMesh::from_triangle_strips(
            square_vertices(),
            &[vec![0, 1, 2], vec![1, 2, 3]],
        )
        .unwrap();

        assert_eq!(mesh.num_triangles(), 2);
    }

    #[test]
    fn triangle_resolves_its_indices_to_vertices() {
        let mesh = TriangleMesh::new(square_vertices(), vec![[0, 1, 3]]).unwrap();

        let t = mesh.triangle(0).expect("the mesh has one face");

        assert_eq!(t.vertices[0], Point3D::from([0.0, 0.0, 0.0]));
        assert_eq!(t.vertices[2], Point3D::from([1.0, 1.0, 0.0]));
    }

    #[test]
    fn triangle_past_the_end_is_none() {
        let mesh = TriangleMesh::new(square_vertices(), vec![[0, 1, 3]]).unwrap();

        assert!(mesh.triangle(1).is_none());
    }

    #[test]
    fn bounds_span_the_vertex_pool() {
        let mesh = TriangleMesh::new(square_vertices(), vec![[0, 1, 2]]).unwrap();

        let (min, max) = mesh.bounds().expect("a non-empty mesh");

        assert_eq!(min, Point3D::from([0.0, 0.0, 0.0]));
        assert_eq!(max, Point3D::from([1.0, 1.0, 0.0]));
    }

    #[test]
    fn from_triangle_makes_a_one_faced_mesh() {
        let t = Triangle3D::new(
            Point3D::from([0.0, 0.0, 0.0]),
            Point3D::from([1.0, 0.0, 0.0]),
            Point3D::from([0.0, 1.0, 0.0]),
        );

        let mesh = TriangleMesh::from_triangle(&t);

        assert_eq!(mesh.num_triangles(), 1);
        assert_eq!(mesh.triangle(0), Some(t));
    }

    #[test]
    fn bounds_of_an_empty_mesh_are_none() {
        assert!(TriangleMesh::default().bounds().is_none());
    }
}
