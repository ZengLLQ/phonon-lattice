use phonon_lattice::dynamics::DynamicalMatrixBuilder;
use phonon_lattice::materials::{GaAsMaterial, GrapheneMaterial, NaClMaterial, SiliconMaterial};
use phonon_lattice::math::{eigensystem_hermitian, Vector3};

#[test]
fn test_silicon_optical_phonon_frequency() {
    let cell = SiliconMaterial::unit_cell();
    let ifc = SiliconMaterial::force_constants();
    let builder = DynamicalMatrixBuilder::new(&cell, &ifc);

    let d_gamma = builder.build_at(Vector3::ZERO);
    let eigen = eigensystem_hermitian(&d_gamma, 1e-13, 100);

    // Silicon optical frequency is ~15.6 THz
    let omega_sq = eigen.eigenvalues[3];
    assert!(omega_sq > 0.0);
    let freq_thz = omega_sq.sqrt() / (2.0 * std::f64::consts::PI * 1e12);

    assert!(
        (freq_thz - 15.6).abs() < 1.0,
        "Silicon optical frequency expected ~15.6 THz, got {:.2} THz",
        freq_thz
    );
}

#[test]
fn test_gaas_loto_splitting() {
    let cell = GaAsMaterial::unit_cell();
    let ifc = GaAsMaterial::force_constants();
    let builder = DynamicalMatrixBuilder::new(&cell, &ifc);

    // Evaluate D(q) approaching Gamma along [1, 0, 0] direction
    let q_small = Vector3::new(1e-4, 0.0, 0.0);
    let d_mat = builder.build_at(q_small);
    let eigen = eigensystem_hermitian(&d_mat, 1e-13, 100);

    // In GaAs: 3 acoustic (~0), 2 TO modes, 1 LO mode higher in frequency
    let freq_to = eigen.eigenvalues[3].sqrt() / (2.0 * std::f64::consts::PI * 1e12);
    let freq_lo = eigen.eigenvalues[5].sqrt() / (2.0 * std::f64::consts::PI * 1e12);

    assert!(
        freq_lo > freq_to,
        "LO mode must exceed TO mode due to polar Born charges: LO={:.2} THz, TO={:.2} THz",
        freq_lo,
        freq_to
    );

    let splitting = freq_lo - freq_to;
    assert!(
        splitting > 0.3,
        "Expected LO-TO splitting in GaAs > 0.3 THz, got {:.3} THz",
        splitting
    );
}

#[test]
fn test_nacl_ionic_loto_splitting() {
    let cell = NaClMaterial::unit_cell();
    let ifc = NaClMaterial::force_constants();
    let builder = DynamicalMatrixBuilder::new(&cell, &ifc);

    let q_small = Vector3::new(0.0, 1e-4, 0.0);
    let d_mat = builder.build_at(q_small);
    let eigen = eigensystem_hermitian(&d_mat, 1e-13, 100);

    let freq_to = eigen.eigenvalues[3].sqrt() / (2.0 * std::f64::consts::PI * 1e12);
    let freq_lo = eigen.eigenvalues[5].sqrt() / (2.0 * std::f64::consts::PI * 1e12);

    assert!(
        freq_lo > freq_to,
        "NaCl ionic LO must exceed TO: LO={:.2} THz, TO={:.2} THz",
        freq_lo,
        freq_to
    );
}

#[test]
fn test_graphene_acoustic_branches() {
    let cell = GrapheneMaterial::unit_cell();
    let ifc = GrapheneMaterial::force_constants();
    let builder = DynamicalMatrixBuilder::new(&cell, &ifc);

    let d_gamma = builder.build_at(Vector3::ZERO);
    let eigen = eigensystem_hermitian(&d_gamma, 1e-13, 100);

    let max_opt = eigen.eigenvalues[5].abs().max(1.0);

    for i in 0..3 {
        let rel_err = eigen.eigenvalues[i].abs() / max_opt;
        assert!(
            rel_err < 1e-10,
            "Graphene acoustic mode #{} at Gamma is non-zero: raw={}, rel_err={}",
            i,
            eigen.eigenvalues[i],
            rel_err
        );
    }
}
