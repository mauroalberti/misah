
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Plane {
    pub point: Point3D,
    pub normal: Vector3D,
}

impl Plane {

    pub fn new(point: Point3D, normal: Vector3D) -> Option<Self> {

        let normal = normal.normalize()?;
        Some(Self { point, normal })
    }

    pub fn from_three_points(
        p0: Point3D,
        p1: Point3D,
        p2: Point3D,
    ) -> Option<Self> {

        let v1 = p0.vector_to(&p1);
        let v2 = p0.vector_to(&p2);

        let normal = v1.cross(&v2).normalize()?;

        Some(Self {
            point: p0,
            normal,
        })
    }

    pub fn signed_distance_to_point(&self, point: &Point3D) -> f64 {

        let v = self.point.vector_to(point);
        self.normal.dot(&v)
    }

    pub fn distance_to_point(&self, point: &Point3D) -> f64 {
        self.signed_distance_to_point(point).abs()
    }

    pub fn contains_point(&self, point: &Point3D, tolerance: f64) -> bool {
        self.distance_to_point(point) <= tolerance
    }

    pub fn coefficients(&self) -> (f64, f64, f64, f64) {

        let [a, b, c] = self.normal.coords;
        let [x0, y0, z0] = self.point.coords;

        let d = -(a * x0 + b * y0 + c * z0);

        (a, b, c, d)
    }
}
