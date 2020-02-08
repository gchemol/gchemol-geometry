// alignment
// :PROPERTIES:
// :header-args: :tangle tests/alignment.rs
// :END:
// ported from: https://theobald.brandeis.edu/QCP/

// [[file:~/Workspace/Programming/gchemol-rs/gchemol-geometry/gchemol-geometry.note::*alignment][alignment:1]]
#[test]
fn test_alignment() {
    use approx::*;
    use gchemol_core::Molecule;
    use gchemol_geometry::Alignment;
    use gchemol_readwrite::prelude::*;
    use vecfx::*;

    // load test molecules
    let mol1 = Molecule::from_file("tests/files/alignment/reference.mol2").expect("alignment reference");
    let mol2 = Molecule::from_file("tests/files/alignment/candidate.mol2").expect("alignment candidate");

    // take the first 5 atoms for superposition
    let reference: Vec<_> = mol1.positions().take(5).collect();
    let candidate: Vec<_> = mol2.positions().take(5).collect();

    // align the candidate onto the reference
    let mut align = Alignment::new(&candidate);
    let sp = align.superpose(&reference, None).unwrap();

    // apply superposition to all atoms
    let new = sp.apply(&candidate);
    assert_relative_eq!(reference.to_matrix(), new.to_matrix(), epsilon = 1e-3);
}
// alignment:1 ends here
