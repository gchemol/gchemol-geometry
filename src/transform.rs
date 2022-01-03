// [[file:../gchemol-geometry.note::7ba19d02][7ba19d02]]
use super::*;
use crate::base::euclidean_distance;

use na::{Rotation3, Vector3};
use nalgebra as na;
// 7ba19d02 ends here

// [[file:../gchemol-geometry.note::fb6473d0][fb6473d0]]
type Points = Vec<Coord3>;

/// Translate all points to a new location
pub fn translate(points: &mut Points, loc: Coord3) {
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
pub fn get_distance_matrix(points: &[Coord3]) -> Vec<Vec<f64>> {
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
pub fn rotate_about_x_axis(points: &Points, angle: f64, center: Coord3) -> Points {
    let axis = Vector3::x_axis();
    let r = Rotation3::from_axis_angle(&axis, angle);

    let mut rpoints = vec![];
    let center = Vector3::from(center);
    for &p in points.iter() {
        let v = Vector3::from(p) - center;
        let t: Coord3 = (r * v + center).into();
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
// fb6473d0 ends here
