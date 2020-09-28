// [[file:../gchemol-geometry.note::*mods][mods:1]]
mod qcprot;
mod quaternion;
// mods:1 ends here

// [[file:../gchemol-geometry.note::*imports][imports:1]]
use gchemol_gut::prelude::*;

use vecfx::*;

type Point3 = [f64; 3];
// imports:1 ends here

// [[file:../gchemol-geometry.note::*base][base:1]]
#[derive(Clone, Copy, Debug)]
pub enum SuperpositionAlgo {
    QCP,
    Quaternion,
}

impl Default for SuperpositionAlgo {
    fn default() -> Self {
        SuperpositionAlgo::QCP
        // SuperpositionAlgo::Quaternion
    }
}

/// The result of alignment defining how to superimpose.
#[derive(Clone, Debug)]
pub struct Superposition {
    /// superpostion rmsd
    pub rmsd: f64,

    /// translation vector
    pub translation: Vector3f,

    /// rotation matrix
    pub rotation_matrix: Matrix3f,
}

impl Superposition {
    /// Apply superposition to other structure `conf`.
    pub fn apply(&self, conf: &[[f64; 3]]) -> Vec<[f64; 3]> {
        let mut res = Vec::with_capacity(conf.len());
        for &v in conf {
            let v = Vector3f::from(v);
            let v = self.rotation_matrix * v + self.translation;
            res.push(v.into());
        }

        // detect NaN floats
        if res.as_flat().iter().any(|x| x.is_nan()) {
            dbg!(&self);
            panic!("found invalid float numbers!");
        }

        res
    }

    /// Apply translation to other structure `conf`.
    pub fn apply_translation(&self, conf: &[Point3]) -> Vec<Point3> {
        let mut res = Vec::with_capacity(conf.len());
        for &v in conf {
            let v = Vector3f::from(v);
            let v = v + self.translation;
            res.push(v.into());
        }

        res
    }

    /// Apply rotation to other structure `conf`.
    pub fn apply_rotation(&self, conf: &[Point3]) -> Vec<Point3> {
        let mut res = Vec::with_capacity(conf.len());
        for &v in conf {
            let v = Vector3f::from(v);
            let v = v + self.rotation_matrix * v;
            res.push(v.into());
        }

        res
    }
}
// base:1 ends here

// [[file:../gchemol-geometry.note::*alignment/deprecated][alignment/deprecated:1]]
#[deprecated(note = "use Superpose instead")]
/// Alignment of candidate structure onto the reference
#[derive(Clone, Debug)]
pub struct Alignment<'a> {
    /// The positions of the candidate structure
    positions: &'a [[f64; 3]],

    /// Select algo
    pub algorithm: SuperpositionAlgo,
}

impl<'a> Alignment<'a> {
    /// Construct from positions of the candidate to be aligned
    pub fn new(positions: &'a [[f64; 3]]) -> Self {
        Alignment {
            positions,
            algorithm: SuperpositionAlgo::default(),
        }
    }

    /// Calculate Root-mean-square deviation of self with the reference coordinates
    ///
    /// Parameters
    /// ----------
    /// * reference: reference coordinates
    /// * weights  : weight of each point
    pub fn rmsd(&self, reference: &[[f64; 3]], weights: Option<&[f64]>) -> Result<f64> {
        // sanity check
        let npts = self.positions.len();
        if reference.len() != npts {
            bail!("points size mismatch!");
        }
        if weights.is_some() && weights.unwrap().len() != npts {
            bail!("weights size mismatch!");
        }

        // calculate rmsd
        let mut ws = 0.0f64;
        for i in 0..npts {
            // take the weight if any, or set it to 1.0
            let wi = weights.map_or_else(|| 1.0, |w| w[i]);
            let dx = wi * (self.positions[i][0] - reference[i][0]);
            let dy = wi * (self.positions[i][1] - reference[i][1]);
            let dz = wi * (self.positions[i][2] - reference[i][2]);

            ws += dx.powi(2) + dy.powi(2) + dz.powi(2);
        }
        let ws = ws.sqrt();

        Ok(ws)
    }

    /// Superpose candidate structure onto reference structure which will be held fixed
    /// Return superposition struct
    ///
    /// Parameters
    /// ----------
    /// * reference: reference coordinates
    /// * weights  : weight of each point
    pub fn superpose(&mut self, reference: &[[f64; 3]], weights: Option<&[f64]>) -> Result<Superposition> {
        // calculate the RMSD & rotational matrix
        let (rmsd, trans, rot) = match self.algorithm {
            SuperpositionAlgo::QCP => self::qcprot::calc_rmsd_rotational_matrix(&reference, &self.positions, weights),
            SuperpositionAlgo::Quaternion => self::quaternion::calc_rmsd_rotational_matrix(&reference, &self.positions, weights),
        };

        // return unit matrix if two structures are already close enough
        let rotation_matrix = if let Some(rot) = rot {
            Matrix3f::from_row_slice(&rot)
        } else {
            Matrix3f::identity()
        };

        // return superimposition result
        let sp = Superposition {
            rmsd,
            translation: trans.into(),
            rotation_matrix,
        };

        Ok(sp)
    }
}
// alignment/deprecated:1 ends here

// [[file:../gchemol-geometry.note::*superimpose][superimpose:1]]
pub type Superimpose<'a> = Superpose<'a>;

/// Superpose of candidate structure onto the reference
#[derive(Clone, Debug)]
pub struct Superpose<'a> {
    /// The positions of the candidate structure
    positions: &'a [[f64; 3]],

    /// Select algo
    pub algorithm: SuperpositionAlgo,
}

impl<'a> Superpose<'a> {
    /// Construct from positions of the candidate to be aligned
    pub fn new(positions: &'a [[f64; 3]]) -> Self {
        Self {
            positions,
            algorithm: SuperpositionAlgo::default(),
        }
    }

    /// Calculate Root-mean-square deviation of self with the reference coordinates
    ///
    /// Parameters
    /// ----------
    /// * reference: reference coordinates
    /// * weights  : weight of each point
    pub fn rmsd(&self, reference: &[[f64; 3]], weights: Option<&[f64]>) -> Result<f64> {
        // sanity check
        let npts = self.positions.len();
        if reference.len() != npts {
            bail!("points size mismatch!");
        }
        if weights.is_some() && weights.unwrap().len() != npts {
            bail!("weights size mismatch!");
        }

        // calculate rmsd
        let mut ws = 0.0f64;
        for i in 0..npts {
            // take the weight if any, or set it to 1.0
            let wi = weights.map_or_else(|| 1.0, |w| w[i]);
            let dx = wi * (self.positions[i][0] - reference[i][0]);
            let dy = wi * (self.positions[i][1] - reference[i][1]);
            let dz = wi * (self.positions[i][2] - reference[i][2]);

            ws += dx.powi(2) + dy.powi(2) + dz.powi(2);
        }
        let ws = ws.sqrt();

        Ok(ws)
    }

    /// Superpose candidate structure onto reference structure which will be held fixed
    /// Return superposition struct
    ///
    /// Parameters
    /// ----------
    /// * reference: reference coordinates
    /// * weights  : weight of each point
    pub fn onto(&mut self, reference: &[[f64; 3]], weights: Option<&[f64]>) -> Superposition {
        // calculate the RMSD & rotational matrix
        let (rmsd, trans, rot) = match self.algorithm {
            SuperpositionAlgo::QCP => self::qcprot::calc_rmsd_rotational_matrix(&reference, &self.positions, weights),
            SuperpositionAlgo::Quaternion => self::quaternion::calc_rmsd_rotational_matrix(&reference, &self.positions, weights),
        };

        // return unit matrix if two structures are already close enough
        let rotation_matrix = if let Some(rot) = rot {
            Matrix3f::from_row_slice(&rot)
        } else {
            Matrix3f::identity()
        };

        // return superimposition result
        Superposition {
            rmsd,
            translation: trans.into(),
            rotation_matrix,
        }
    }
}
// superimpose:1 ends here

// [[file:../gchemol-geometry.note::*test][test:1]]
#[test]
fn test_alignment() {
    use approx::*;

    // fragment a
    let (reference, candidate, weights) = qcprot::prepare_test_data();

    // construct alignment for superimposition
    let sp = Superpose::new(&candidate).onto(&reference, Some(&weights));
    let rot = sp.rotation_matrix;

    // validation
    let rot_expected = Matrix3f::from_row_slice(&[
        0.77227551,
        0.63510272,
        -0.01533190,
        -0.44544846,
        0.52413614,
        -0.72584914,
        -0.45295276,
        0.56738509,
        0.68768304,
    ]);
    assert_relative_eq!(rot_expected, rot, epsilon = 1e-4);
}
// test:1 ends here

// [[file:../gchemol-geometry.note::*test][test:2]]
#[test]
fn test_alignment_hcn() {
    use vecfx::*;

    let positions_ref = [
        [0.83334699, 0.716865, 0.0],
        [-0.5486581, -0.35588, 0.0],
        [-0.2855828, 1.036928, 0.0],
    ];

    let positions_can = [[-0.634504, -0.199638, -0.0], [0.970676, 0.670662, 0.0], [-0.337065, 0.926883, 0.0]];

    let weights = vec![0.0001; 3];
    let sp = Superpose::new(&positions_can).onto(&positions_ref, Some(&weights));
    assert_relative_eq!(sp.rmsd, 0.0614615, epsilon = 1e-4);

    let t = Vector3f::from([0.423160235, 0.2715202, 0.0]);
    let r = Matrix3f::from_column_slice(&[-0.4167190, -0.9090353, 0.0, -0.909035, 0.4167190, 0.0, 0.0, 0.0, -1.0]);
    assert_relative_eq!(sp.rotation_matrix, r, epsilon = 1e-4);
}
// test:2 ends here
