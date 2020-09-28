// [[file:../gchemol-geometry.note::*docs][docs:1]]
//! # Example
//!
//! ```no_run
//! use gchemol_geometry::Superpose;
//! use gchemol_core::Molecule;
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

// [[file:../gchemol-geometry.note::*mods][mods:1]]
mod alignment;
mod base;
mod random;
mod transform;
// mods:1 ends here

// [[file:../gchemol-geometry.note::*pub][pub:1]]
pub use crate::alignment::*;

#[cfg(feature="adhoc")]
pub use crate::base::*;

#[cfg(feature="adhoc")]
pub use crate::random::*;

#[cfg(feature="adhoc")]
pub use crate::transform::*;

#[cfg(test)]
#[macro_use]
extern crate approx;
// pub:1 ends here
