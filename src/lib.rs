// [[file:../gchemol-geometry.note::*docs][docs:1]]
//! # Example
//!
//! ```no_run
//! use gchemol::geom::Superpose;
//! use gchemol::Molecule;
//! use gchemol::prelude::*;
//! 
//! // load test molecules
//! let mol1 = Molecule::from_file("tests/files/alignment/reference.mol2").unwrap();
//! let mol2 = Molecule::from_file("tests/files/alignment/candidate.mol2").unwrap();
//! 
//! // take the first 5 atoms for superposition
//! let reference: Vec<_> = mol1.positions().take(5).collect();
//! let candidate: Vec<_> = mol2.positions().take(5).collect();
//! 
//! // align the candidate onto the reference
//! let sp = Superpose::new(&candidate).onto(&reference, None);
//! 
//! // apply superposition to all atoms
//! let superimposed_structure = sp.apply(&candidate);
//! 
//! // apply translation only
//! let translated_structure = sp.apply_translation(&candidate);
//! ```
// docs:1 ends here

// [[file:../gchemol-geometry.note::f065136b][f065136b]]
use gchemol_gut::prelude::*;
use vecfx::*;
// f065136b ends here

// [[file:../gchemol-geometry.note::a70e28c8][a70e28c8]]
mod alignment;
mod base;
mod random;
mod traits;
mod transform;
// a70e28c8 ends here

// [[file:../gchemol-geometry.note::62451bc9][62451bc9]]
/// Three-dimensional Cartesian coordinates
pub type Coord3 = [f64; 3];

pub use crate::alignment::*;
pub use crate::base::*;

#[cfg(feature = "adhoc")]
pub use crate::random::*;

#[cfg(feature = "adhoc")]
pub use crate::transform::*;

pub mod prelude {
    pub use crate::traits::*;
}
// 62451bc9 ends here
