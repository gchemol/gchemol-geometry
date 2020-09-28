// [[file:../gchemol-geometry.note::*imports][imports:1]]
use gchemol_gut::prelude::*;
use vecfx::*;
// imports:1 ends here

// [[file:../gchemol-geometry.note::*functions][functions:1]]
#[inline]
/// Return Cartesian distance between two points in 3D space.
pub fn euclidean_distance(p1: [f64; 3], p2: [f64; 3]) -> f64 {
    let mut d2 = 0.0;
    for v in 0..3 {
        let dv = p2[v] - p1[v];
        d2 += dv * dv;
    }

    d2.sqrt()
}

/// Return the geometric center
pub fn weighted_center_of_geometry(positions: &[[f64; 3]], weights: &[f64]) -> Result<Vector3f> {
    let npts = positions.len();
    let mut pc = [0.0; 3];

    // sanity check
    if npts != weights.len() {
        bail!("size inconsistent!");
    }

    // deviding by zero?
    let wsum = weights.sum();
    if wsum < 1e-6 {
        error!("weird weight sum: {:?}", wsum);
    }

    for i in 0..npts {
        for j in 0..3 {
            pc[j] += weights[i] * positions[i][j];
        }
    }

    for i in 0..3 {
        pc[i] /= wsum;
    }

    Ok(Vector3f::from(pc))
}
// functions:1 ends here

// [[file:../gchemol-geometry.note::*test][test:1]]
#[test]
fn test_weighted_center_of_geometry() {
    use approx::*;

    // points
    let frag = vec![
        [-2.803, -15.373, 24.556],
        [0.893, -16.062, 25.147],
        [1.368, -12.371, 25.885],
        [-1.651, -12.153, 28.177],
        [-0.440, -15.218, 30.068],
        [2.551, -13.273, 31.372],
        [0.105, -11.330, 33.567],
    ];

    // weights
    let natoms = frag.len();
    let masses: Vec<_> = (0..natoms).map(|v| v as f64 + 1.0).collect();

    // expected results
    let expected = Vector3f::new(0.3687142857142857, -13.15214285714286, 29.955499999999997);
    let pc = weighted_center_of_geometry(&frag, &masses).expect("geometry: com");
    assert_relative_eq!(pc, expected, epsilon = 1e-6);
}
// test:1 ends here
