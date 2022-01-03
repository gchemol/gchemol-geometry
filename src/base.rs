// [[file:../gchemol-geometry.note::ee1c9917][ee1c9917]]
use super::*;
// ee1c9917 ends here

// [[file:../gchemol-geometry.note::3dde6602][3dde6602]]
#[inline]
/// Return Cartesian distance between two points in 3D space.
pub fn euclidean_distance(p1: Coord3, p2: Coord3) -> f64 {
    let mut d2 = 0.0;
    for v in 0..3 {
        let dv = p2[v] - p1[v];
        d2 += dv * dv;
    }

    d2.sqrt()
}

// FIXME: when sum of weight is too large
/// Return the weighted geometric center
pub fn weighted_center_of_geometry(positions: &[Coord3], weights: &[f64]) -> Coord3 {
    let npts = positions.len();
    assert_eq!(npts, weights.len(), "array size mismatch between positions and weights");

    // deviding by zero?
    let wsum = weights.sum();
    assert_ne!(wsum, 0.0, "invalid sum of weights");

    let mut pc = [0.0; 3];
    for i in 0..npts {
        for j in 0..3 {
            pc[j] += weights[i] * positions[i][j];
        }
    }

    for i in 0..3 {
        pc[i] /= wsum;
    }

    pc
}
// 3dde6602 ends here

// [[file:../gchemol-geometry.note::30f52d4f][30f52d4f]]
#[test]
fn test_weighted_center_of_geometry() {
    use vecfx::approx::assert_relative_eq;

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
    let pc = weighted_center_of_geometry(&frag, &masses);
    assert_relative_eq!(Vector3f::from(pc), expected, epsilon = 1e-6);
}
// 30f52d4f ends here
