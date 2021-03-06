// imports

// [[file:~/Workspace/Programming/gchemol-rs/gchemol-geometry/gchemol-geometry.note::*imports][imports:1]]
use crate::base::euclidean_distance;

use vecfx::nalgebra as na;
use vecfx::*;

use na::{Rotation3, Vector3};
// imports:1 ends here

// impl

// [[file:~/Workspace/Programming/gchemol-rs/gchemol-geometry/gchemol-geometry.note::*impl][impl:1]]
type Point3 = [f64; 3];
type Points = Vec<Point3>;

/// Translate all points to a new location
pub fn translate(points: &mut Points, loc: Point3) {
    for i in 0..points.len() {
        for v in 0..3 {
            points[i][v] += loc[v];
        }
    }
}

/// check if any pair of points come too close
pub fn close_contact(points: &Points) -> bool {
    let cutoff = 0.4;

    let npts = points.len();
    for i in 0..npts {
        for j in (i + 1)..npts {
            let p1 = points[i];
            let p2 = points[j];
            let dx = p2[0] - p1[0];
            let dy = p2[1] - p1[1];
            let dz = p2[2] - p1[2];
            let d2 = dx * dx + dy * dy + dz * dz;
            if d2 <= cutoff {
                return true;
            }
        }
    }

    false
}

/// Return all distances between any pair of points
pub fn get_distance_matrix(points: &[Point3]) -> Vec<Vec<f64>> {
    let npts = points.len();

    // fill distance matrix
    let mut distmat = vec![];
    for i in 0..npts {
        let mut dijs = vec![];
        for j in 0..npts {
            let dij = euclidean_distance(points[i], points[j]);
            dijs.push(dij);
        }
        distmat.push(dijs);
    }

    distmat
}

/// rotate coordinates about x axis in radian
pub fn rotate_about_x_axis(points: &Points, angle: f64, center: Point3) -> Points {
    let axis = Vector3::x_axis();
    let r = Rotation3::from_axis_angle(&axis, angle);

    let mut rpoints = vec![];
    let center = Vector3::from(center);
    for &p in points.iter() {
        let v = Vector3::from(p) - center;
        let t: Point3 = (r * v + center).into();
        rpoints.push(t);
    }

    rpoints
}

/// Return mirror inverted structure
pub fn mirror_invert(positions: &[[f64; 3]]) -> Vector3fVec {
    let m = positions.to_matrix();
    let r = na::Matrix3::from_diagonal(&[1.0, 1.0, -1.0].into());
    r * m
}

/// Return point inverted structure
pub fn point_invert(positions: &[[f64; 3]]) -> Vector3fVec {
    let m = positions.to_matrix();
    let r = na::Matrix3::from_diagonal(&[-1.0, -1.0, -1.0].into());
    r * m
}
// impl:1 ends here
