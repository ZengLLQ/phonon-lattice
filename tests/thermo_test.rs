use phonon_lattice::crystal::ReciprocalLattice;
use phonon_lattice::dynamics::DynamicalMatrixBuilder;
use phonon_lattice::materials::SiliconMaterial;
use phonon_lattice::math::eigensystem_hermitian;
use phonon_lattice::thermodynamics::ThermodynamicCalculator;

#[test]
fn test_dulong_petit_classical_limit() {
    let cell = SiliconMaterial::unit_cell();
    let ifc = SiliconMaterial::force_constants();
    let recip = ReciprocalLattice::from_unit_cell(&cell);
    let builder = DynamicalMatrixBuilder::new(&cell, &ifc);

    let q_pts = recip.generate_monkhorst_pack_grid(4, 4, 4);
    let n_q = q_pts.len();
    let num_branches = 3 * cell.num_atoms();

    let mut freqs = Vec::with_capacity(n_q * num_branches);
    let mut weights = Vec::with_capacity(n_q * num_branches);
    let w_single = 1.0 / n_q as f64;

    for q in &q_pts {
        let cart_q = recip.fractional_to_cartesian(*q);
        let d = builder.build_at(cart_q);
        let eig = eigensystem_hermitian(&d, 1e-12, 100);
        for &val in &eig.eigenvalues {
            let omega = if val > 0.0 { val.sqrt() } else { 0.0 };
            freqs.push(omega);
            weights.push(w_single);
        }
    }

    let thermo =
        ThermodynamicCalculator::new(freqs, weights, cell.num_atoms(), cell.total_mass_kg());

    // At high T (5000 K), Cv must approach 3 * N_atoms * R = 3 * 2 * 8.314 = 49.887 J/(mol K)
    let state_high_t = thermo.evaluate_at_temperature(5000.0);
    let expected_dulong_petit = 3.0 * (cell.num_atoms() as f64) * 8.314462618;

    assert!(
        (state_high_t.heat_capacity_cv_j_per_mol_k - expected_dulong_petit).abs() < 1.0,
        "Dulong-Petit limit failed: expected {:.2}, got {:.2}",
        expected_dulong_petit,
        state_high_t.heat_capacity_cv_j_per_mol_k
    );
}

#[test]
fn test_thermodynamic_monotonicity_and_zpe() {
    let cell = SiliconMaterial::unit_cell();
    let ifc = SiliconMaterial::force_constants();
    let recip = ReciprocalLattice::from_unit_cell(&cell);
    let builder = DynamicalMatrixBuilder::new(&cell, &ifc);

    let q_pts = recip.generate_monkhorst_pack_grid(4, 4, 4);
    let n_q = q_pts.len();
    let num_branches = 3 * cell.num_atoms();

    let mut freqs = Vec::with_capacity(n_q * num_branches);
    let mut weights = Vec::with_capacity(n_q * num_branches);
    let w_single = 1.0 / n_q as f64;

    for q in &q_pts {
        let cart_q = recip.fractional_to_cartesian(*q);
        let d = builder.build_at(cart_q);
        let eig = eigensystem_hermitian(&d, 1e-12, 100);
        for &val in &eig.eigenvalues {
            let omega = if val > 0.0 { val.sqrt() } else { 0.0 };
            freqs.push(omega);
            weights.push(w_single);
        }
    }

    let thermo =
        ThermodynamicCalculator::new(freqs, weights, cell.num_atoms(), cell.total_mass_kg());

    let state_low = thermo.evaluate_at_temperature(10.0);
    let state_mid = thermo.evaluate_at_temperature(300.0);
    let state_high = thermo.evaluate_at_temperature(800.0);

    // Heat capacity monotonic increase
    assert!(state_low.heat_capacity_cv_j_per_mol_k < state_mid.heat_capacity_cv_j_per_mol_k);
    assert!(state_mid.heat_capacity_cv_j_per_mol_k < state_high.heat_capacity_cv_j_per_mol_k);

    // ZPE must be strictly positive
    assert!(state_mid.zero_point_energy_kj_per_mol > 1.0);
}
