use phonon_lattice::dynamics::DynamicalMatrixBuilder;
use phonon_lattice::materials::SiliconMaterial;
use phonon_lattice::math::{eigensystem_hermitian, Vector3};

#[test]
fn test_acoustic_sum_rule_at_gamma() {
    let cell = SiliconMaterial::unit_cell();
    let ifc = SiliconMaterial::force_constants();
    let builder = DynamicalMatrixBuilder::new(&cell, &ifc);

    // At Gamma (q = 0, 0, 0), the 3 lowest acoustic modes must have eigenvalue == 0
    let d_gamma = builder.build_at(Vector3::ZERO);
    let eigen = eigensystem_hermitian(&d_gamma, 1e-14, 100);

    assert_eq!(eigen.eigenvalues.len(), 6);

    let max_opt = eigen.eigenvalues[5].abs().max(1.0);

    // First 3 eigenvalues must be vanishingly small compared to optical modes (relative error < 1e-12)
    for i in 0..3 {
        let rel_err = eigen.eigenvalues[i].abs() / max_opt;
        assert!(
            rel_err < 1e-10,
            "Acoustic mode #{} at Gamma is non-zero: raw={}, rel_err={}",
            i,
            eigen.eigenvalues[i],
            rel_err
        );
    }

    // Remaining 3 optical modes in diamond silicon are triply degenerate at Gamma
    let opt0 = eigen.eigenvalues[3];
    let opt1 = eigen.eigenvalues[4];
    let opt2 = eigen.eigenvalues[5];

    assert!(
        (opt0 - opt1).abs() < 1e-6,
        "Optical modes split: {} vs {}",
        opt0,
        opt1
    );
    assert!(
        (opt1 - opt2).abs() < 1e-6,
        "Optical modes split: {} vs {}",
        opt1,
        opt2
    );
}

#[test]
fn test_dynamical_matrix_hermiticity_at_arbitrary_q() {
    let cell = SiliconMaterial::unit_cell();
    let ifc = SiliconMaterial::force_constants();
    let builder = DynamicalMatrixBuilder::new(&cell, &ifc);

    let q = Vector3::new(0.35, 0.12, 0.78);
    let d = builder.build_at(q);

    let dim = d.len();
    for i in 0..dim {
        for j in 0..dim {
            let diff = d[i][j] - d[j][i].conj();
            assert!(
                diff.norm() < 1e-12,
                "Hermiticity violated at ({}, {}): d[i][j]={:?}, d[j][i]={:?}",
                i,
                j,
                d[i][j],
                d[j][i]
            );
        }
    }
}
