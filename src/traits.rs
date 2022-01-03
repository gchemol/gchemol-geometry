// [[file:../gchemol-geometry.note::c1257b7a][c1257b7a]]
use super::*;
// c1257b7a ends here

// [[file:../gchemol-geometry.note::23e4530d][23e4530d]]
/// A trait provides useful methods for [f64; 3] type.
pub trait GeometryCoord3Ext {
    /// return the distance to other point
    fn distance(&self, other: Self) -> f64;

    /// return the angle between the tree vector points: va(self), vb, vc
    fn angle(&self, pb: Self, pc: Self) -> f64;

    /// return the torsion angle between the four vector points: va, vb, vc, vd
    fn torsion(&self, pb: Self, pc: Self, pd: Self) -> f64;
}

impl GeometryCoord3Ext for Coord3 {
    /// return the distance to other point
    fn distance(&self, other: Self) -> f64 {
        euclidean_distance(*self, other.into())
    }

    /// return the angle between three points: p0(self), p1, p2
    fn angle(&self, p1: Self, p2: Self) -> f64 {
        let a = self.as_vector_slice();
        let b = p1.as_vector_slice();
        let c = p2.as_vector_slice();
        let ab = a - b;
        let cb = c - b;
        ab.angle(&cb)
    }

    /// return the torsion angle from four points: p0(self), p1, p2, p3
    ///
    /// The dihedral measures the rotation around p1-p2
    ///
    ///    p0---->p1
    ///           \
    ///            \
    ///             p2---->p3
    ///
    /// The dihedral angle is restricted to the range -π <= x <= π.
    ///
    /// Reference
    /// ---------
    /// http://stackoverflow.com/a/34245697
    fn torsion(&self, p1: Self, p2: Self, p3: Self) -> f64 {
        let p0 = self.as_vector_slice();
        let p1 = p1.as_vector_slice();
        let p2 = p2.as_vector_slice();
        let p3 = p3.as_vector_slice();

        let b0 = p0 - p1;
        let b1 = (p2 - p1).normalize();
        let b2 = p3 - p2;

        // v = projection of b0 onto plane perpendicular to b1
        // w = projection of b2 onto plane perpendicular to b1
        // let v = b0.vector_rejection(&b1);
        // let w = b2.vector_rejection(&b1);
        let v = &b0 - b0.dot(&b1) * &b1;
        let w = &b2 - b2.dot(&b1) * &b1;

        // angle between v and w in a plane is the torsion angle
        // v and w may not be normalized but that's fine since tan is y/x
        let x = v.dot(&w);
        let y = w.dot(&b1.cross(&v));
        y.atan2(x)
    }
}
// 23e4530d ends here

// [[file:../gchemol-geometry.note::86305981][86305981]]
/// A trait provides useful methods for slice of [f64; 3] type.
pub trait GeometryCoord3SliceExt {
    /// Return the center of geometry (COG)
    fn center_of_geometry(&self) -> Coord3;

    /// Return the center of geometry
    ///
    /// https://en.wikipedia.org/wiki/Centroid
    fn centroid(&self) -> Coord3 {
        self.center_of_geometry()
    }

    /// Return weighted center of geometry (COM)
    fn center_of_mass(&self, masses: &[f64]) -> Coord3;

    /// Apply mirror inversion.
    fn mirror_invert(&mut self);

    /// Apply point inversion.
    fn point_invert(&mut self);
}

impl GeometryCoord3SliceExt for [Coord3] {
    /// Return weighted center of geometry (COM)
    fn center_of_mass(&self, masses: &[f64]) -> Coord3 {
        weighted_center_of_geometry(&self, &masses)
    }

    /// Return the center of geometry (COG)
    fn center_of_geometry(&self) -> Coord3 {
        let n = self.len();
        let weights: Vec<_> = (0..n).map(|_| 1.0).collect();
        weighted_center_of_geometry(&self, &weights)
    }

    /// Apply mirror inversion along z-axis.
    fn mirror_invert(&mut self) {
        for p in self.iter_mut() {
            p[2] *= -1.0;
        }
    }

    /// Apply mirror inversion around the center.
    fn point_invert(&mut self) {
        for p in self.iter_mut() {
            p[0] *= -1.0;
            p[1] *= -1.0;
            p[2] *= -1.0;
        }
    }
}
// 86305981 ends here

// [[file:../gchemol-geometry.note::fb603613][fb603613]]
#[test]
fn test_point3_math() {
    use approx::*;
    use std::f64::consts::PI;

    let p1 = [-2.0291693, 1.7660114, 0.0000000];
    let p2 = [-1.6725149, 0.7572014, 0.0000000];
    let p3 = [-1.6724965, 2.2704096, 0.8736515];
    let p4 = [-1.6724965, 2.2704096, -0.8736515];
    let p5 = [-3.0991693, 1.7660246, 0.0000000];

    assert_relative_eq!(p1.distance(p2), 1.07, epsilon = 1e-4);
    assert_relative_eq!(p1.angle(p2, p3).to_degrees(), 35.264, epsilon = 1e-3);
    assert_relative_eq!(p2.torsion(p3, p1, p5).to_degrees(), 120.00, epsilon = 1e-1);
    assert_relative_eq!(p2.torsion(p3, p4, p5).to_degrees(), 70.529, epsilon = 1e-1);
    assert_relative_eq!(p1.torsion(p2, p3, p4).to_degrees(), -35.246, epsilon = 1e-1);
}
// fb603613 ends here
