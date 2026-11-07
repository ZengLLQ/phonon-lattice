use phonon_lattice::crystal::{Atom, ReciprocalLattice, UnitCell};
use phonon_lattice::math::Vector3;
use std::f64::consts::PI;

#[test]
fn test_cubic_unit_cell_volume_and_metric() {
    let a = 5.0; // Angstroms
    let a1 = Vector3::new(a, 0.0, 0.0);
    let a2 = Vector3::new(0.0, a, 0.0);
    let a3 = Vector3::new(0.0, 0.0, a);
    let atoms = vec![Atom::new("Si", 28.0855, Vector3::ZERO)];

    let cell = UnitCell::new(a1, a2, a3, atoms);
    assert!((cell.volume() - 125.0).abs() < 1e-12);
    assert!((cell.volume_m3() - 125.0e-30).abs() < 1e-40);
}

#[test]
fn test_reciprocal_lattice_orthogonality() {
    // FCC primitive cell
    let a = 5.431 * 0.5;
    let a1 = Vector3::new(0.0, a, a);
    let a2 = Vector3::new(a, 0.0, a);
    let a3 = Vector3::new(a, a, 0.0);
    let atoms = vec![Atom::new("Si", 28.0855, Vector3::ZERO)];
    let cell = UnitCell::new(a1, a2, a3, atoms);

    let recip = ReciprocalLattice::from_unit_cell(&cell);

    // Assert a_i . b_j = 2 * pi * delta_ij
    let direct = [cell.a1, cell.a2, cell.a3];
    let reciprocal = [recip.b1, recip.b2, recip.b3];

    for i in 0..3 {
        for j in 0..3 {
            let dot = direct[i].dot(&reciprocal[j]);
            let expected = if i == j { 2.0 * PI } else { 0.0 };
            assert!(
                (dot - expected).abs() < 1e-10,
                "Orthogonality failed for ({}, {}): got {}, expected {}",
                i,
                j,
                dot,
                expected
            );
        }
    }
}

#[test]
fn test_fractional_cartesian_roundtrip() {
    let a = 5.0;
    let cell = UnitCell::new(
        Vector3::new(a, 0.0, 0.0),
        Vector3::new(0.0, a, 0.0),
        Vector3::new(0.0, 0.0, a),
        vec![],
    );

    let frac = Vector3::new(0.35, 0.72, 0.18);
    let cart = cell.fractional_to_cartesian(frac);
    let back_frac = cell.cartesian_to_fractional(cart).expect("Inverse exists");

    assert!((back_frac.x - frac.x).abs() < 1e-12);
    assert!((back_frac.y - frac.y).abs() < 1e-12);
    assert!((back_frac.z - frac.z).abs() < 1e-12);
}
