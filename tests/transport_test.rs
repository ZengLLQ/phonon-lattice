use phonon_lattice::crystal::ReciprocalLattice;
use phonon_lattice::dynamics::DynamicalMatrixBuilder;
use phonon_lattice::materials::SiliconMaterial;
use phonon_lattice::thermodynamics::DebyeModel;
use phonon_lattice::transport::{PhononScatteringModel, ThermalConductivitySolver};

#[test]
fn test_cubic_silicon_thermal_conductivity_isotropy() {
    let cell = SiliconMaterial::unit_cell();
    let ifc = SiliconMaterial::force_constants();
    let recip = ReciprocalLattice::from_unit_cell(&cell);
    let builder = DynamicalMatrixBuilder::new(&cell, &ifc);

    let debye = DebyeModel::compute(8430.0, 5840.0, cell.num_atoms(), cell.volume_m3());

    let scattering = PhononScatteringModel::new(
        1e-4,
        2.01e-4,
        1.10,
        debye.mean_sound_velocity_m_s,
        cell.total_mass_kg() / cell.num_atoms() as f64,
        debye.debye_temperature_k,
        debye.debye_frequency_thz * 2.0 * std::f64::consts::PI * 1e12,
    );

    let solver = ThermalConductivitySolver::new(&builder, &recip, scattering);
    let res = solver.solve_at_temperature(300.0, (4, 4, 4));

    let kxx = res.conductivity_tensor_w_per_m_k.data[0][0];
    let kyy = res.conductivity_tensor_w_per_m_k.data[1][1];
    let kzz = res.conductivity_tensor_w_per_m_k.data[2][2];

    // Must be strictly positive
    assert!(kxx > 10.0, "kxx must be positive: got {}", kxx);
    assert!(kyy > 10.0, "kyy must be positive: got {}", kyy);
    assert!(kzz > 10.0, "kzz must be positive: got {}", kzz);

    // Diagonal elements in cubic lattice must be approximately equal
    let diff_xy = (kxx - kyy).abs() / kxx;
    let diff_xz = (kxx - kzz).abs() / kxx;
    assert!(
        diff_xy < 0.25,
        "Anisotropy between xx and yy too large: {:.3}",
        diff_xy
    );
    assert!(
        diff_xz < 0.25,
        "Anisotropy between xx and zz too large: {:.3}",
        diff_xz
    );

    // Off-diagonal elements should be much smaller than diagonals
    let kxy = res.conductivity_tensor_w_per_m_k.data[0][1];
    assert!(
        kxy.abs() < kxx * 0.15,
        "Off-diagonal kxy too large: {}",
        kxy
    );
}

#[test]
fn test_thermal_conductivity_umklapp_temperature_decay() {
    let cell = SiliconMaterial::unit_cell();
    let ifc = SiliconMaterial::force_constants();
    let recip = ReciprocalLattice::from_unit_cell(&cell);
    let builder = DynamicalMatrixBuilder::new(&cell, &ifc);

    let debye = DebyeModel::compute(8430.0, 5840.0, cell.num_atoms(), cell.volume_m3());

    let scattering = PhononScatteringModel::new(
        1e-4,
        2.01e-4,
        1.10,
        debye.mean_sound_velocity_m_s,
        cell.total_mass_kg() / cell.num_atoms() as f64,
        debye.debye_temperature_k,
        debye.debye_frequency_thz * 2.0 * std::f64::consts::PI * 1e12,
    );

    let solver = ThermalConductivitySolver::new(&builder, &recip, scattering);

    let res_300 = solver.solve_at_temperature(300.0, (4, 4, 4));
    let res_600 = solver.solve_at_temperature(600.0, (4, 4, 4));

    // High-temperature Umklapp scattering increases ~ T, leading to kappa ~ 1/T decay
    assert!(
        res_300.isotropic_conductivity_w_per_m_k > res_600.isotropic_conductivity_w_per_m_k,
        "Thermal conductivity must decrease from 300K ({:.2}) to 600K ({:.2})",
        res_300.isotropic_conductivity_w_per_m_k,
        res_600.isotropic_conductivity_w_per_m_k
    );
}
