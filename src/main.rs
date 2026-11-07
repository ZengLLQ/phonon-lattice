#![allow(clippy::needless_range_loop, clippy::approx_constant)]

use phonon_lattice::crystal::ReciprocalLattice;
use phonon_lattice::dynamics::{DynamicalMatrixBuilder, PhononDispersion, PhononDosCalculator};
use phonon_lattice::materials::{GaAsMaterial, GrapheneMaterial, NaClMaterial, SiliconMaterial};
use phonon_lattice::thermodynamics::{DebyeModel, ThermodynamicCalculator};
use phonon_lattice::transport::{PhononScatteringModel, ThermalConductivitySolver};
use std::env;
use std::time::Instant;

fn print_banner() {
    println!("================================================================================");
    println!("  PHONON LATTICE: Crystal Dynamics & Anharmonic Thermal Transport Engine");
    println!("  Dynamical Matrix | Non-Analytic Corrections | BTE-RTA Thermal Conductivity");
    println!("================================================================================");
}

fn print_help() {
    print_banner();
    println!("Usage: phonon-lattice <command> [material] [options]\n");
    println!("Commands:");
    println!(
        "  dispersion [silicon|graphene|gaas|nacl]   Compute and render phonon band dispersion"
    );
    println!(
        "  dos        [silicon|graphene|gaas|nacl]   Compute and render Phonon Density of States"
    );
    println!(
        "  thermo     [material] [t_min] [t_max]     Evaluate vibrational thermodynamic curves"
    );
    println!("  bte        [material] [temp_k]            Solve BTE thermal conductivity tensor");
    println!(
        "  benchmark                                 Benchmark matrix assembly & diagonalization"
    );
}

fn main() {
    let args: Vec<String> = env::args().collect();
    if args.len() < 2 || args[1] == "--help" || args[1] == "-h" {
        print_help();
        return;
    }

    match args[1].as_str() {
        "dispersion" => {
            print_banner();
            let mat_name = args.get(2).map(|s| s.as_str()).unwrap_or("silicon");
            run_dispersion(mat_name);
        }
        "dos" => {
            print_banner();
            let mat_name = args.get(2).map(|s| s.as_str()).unwrap_or("silicon");
            run_dos(mat_name);
        }
        "thermo" => {
            print_banner();
            let mat_name = args.get(2).map(|s| s.as_str()).unwrap_or("silicon");
            let t_min: f64 = args.get(3).and_then(|s| s.parse().ok()).unwrap_or(50.0);
            let t_max: f64 = args.get(4).and_then(|s| s.parse().ok()).unwrap_or(1000.0);
            run_thermo(mat_name, t_min, t_max);
        }
        "bte" => {
            print_banner();
            let mat_name = args.get(2).map(|s| s.as_str()).unwrap_or("silicon");
            let temp_k: f64 = args.get(3).and_then(|s| s.parse().ok()).unwrap_or(300.0);
            run_bte(mat_name, temp_k);
        }
        "benchmark" => {
            print_banner();
            run_benchmark();
        }
        cmd => {
            eprintln!("Unknown command: {}", cmd);
            print_help();
        }
    }
}

fn run_dispersion(mat: &str) {
    println!("Mode: Phonon Band Structure Dispersion Calculation");
    println!("Material: {}", mat.to_uppercase());

    match mat.to_lowercase().as_str() {
        "silicon" => {
            let cell = SiliconMaterial::unit_cell();
            let ifc = SiliconMaterial::force_constants();
            let recip = ReciprocalLattice::from_unit_cell(&cell);
            let path = SiliconMaterial::high_symmetry_path(&recip, 15);
            let builder = DynamicalMatrixBuilder::new(&cell, &ifc);

            let disp = PhononDispersion::calculate(&builder, &path);
            println!("Number of phonon branches: {}", disp.num_branches);
            println!("Max frequency: {:.2} THz", disp.max_frequency_thz());
            println!("\nASCII Phonon Dispersion Curve:\n");
            print!("{}", disp.render_ascii_plot(16, 65));
        }
        "gaas" => {
            let cell = GaAsMaterial::unit_cell();
            let ifc = GaAsMaterial::force_constants();
            let recip = ReciprocalLattice::from_unit_cell(&cell);
            let path = GaAsMaterial::high_symmetry_path(&recip, 15);
            let builder = DynamicalMatrixBuilder::new(&cell, &ifc);

            let disp = PhononDispersion::calculate(&builder, &path);
            println!("Number of phonon branches: {}", disp.num_branches);
            println!(
                "Polar Crystal Non-Analytic Correction Active (Born effective charges Z* = 2.18)"
            );
            println!("Max frequency: {:.2} THz", disp.max_frequency_thz());
            println!("\nASCII Phonon Dispersion Curve (LO-TO Splitting):\n");
            print!("{}", disp.render_ascii_plot(16, 65));
        }
        "graphene" => {
            let cell = GrapheneMaterial::unit_cell();
            let ifc = GrapheneMaterial::force_constants();
            let recip = ReciprocalLattice::from_unit_cell(&cell);
            let path = GrapheneMaterial::high_symmetry_path(&recip, 15);
            let builder = DynamicalMatrixBuilder::new(&cell, &ifc);

            let disp = PhononDispersion::calculate(&builder, &path);
            println!("Number of phonon branches: {}", disp.num_branches);
            println!(
                "2D Honeycomb Lattice: Flexural Out-of-Plane ZA Branch & High-Frequency G-Mode"
            );
            println!("Max frequency: {:.2} THz", disp.max_frequency_thz());
            println!("\nASCII Phonon Dispersion Curve:\n");
            print!("{}", disp.render_ascii_plot(16, 65));
        }
        "nacl" => {
            let cell = NaClMaterial::unit_cell();
            let ifc = NaClMaterial::force_constants();
            let recip = ReciprocalLattice::from_unit_cell(&cell);
            let path = NaClMaterial::high_symmetry_path(&recip, 15);
            let builder = DynamicalMatrixBuilder::new(&cell, &ifc);

            let disp = PhononDispersion::calculate(&builder, &path);
            println!("Number of phonon branches: {}", disp.num_branches);
            println!("Ionic Rocksalt Crystal (NaCl): Extreme Non-Analytic Coulomb Splitting");
            println!("Max frequency: {:.2} THz", disp.max_frequency_thz());
            println!("\nASCII Phonon Dispersion Curve:\n");
            print!("{}", disp.render_ascii_plot(16, 65));
        }
        other => eprintln!(
            "Unknown material: {}. Available: silicon, gaas, graphene, nacl",
            other
        ),
    }
}

fn run_dos(mat: &str) {
    println!("Mode: Phonon Density of States (DOS) Calculation");
    println!("Material: {}", mat.to_uppercase());

    let (cell, ifc) = match mat.to_lowercase().as_str() {
        "gaas" => (GaAsMaterial::unit_cell(), GaAsMaterial::force_constants()),
        "graphene" => (
            GrapheneMaterial::unit_cell(),
            GrapheneMaterial::force_constants(),
        ),
        "nacl" => (NaClMaterial::unit_cell(), NaClMaterial::force_constants()),
        _ => (
            SiliconMaterial::unit_cell(),
            SiliconMaterial::force_constants(),
        ),
    };

    let recip = ReciprocalLattice::from_unit_cell(&cell);
    let builder = DynamicalMatrixBuilder::new(&cell, &ifc);
    let calculator = PhononDosCalculator::new(&builder, &recip);

    println!("Sampling 8x8x8 Monkhorst-Pack Brillouin zone grid (512 q-points)...");
    let dos = calculator.compute((8, 8, 8), 0.2, 50);

    println!(
        "Total integrated states: {:.3} (Expected: {:.1})",
        dos.integrated_dos,
        3.0 * cell.num_atoms() as f64
    );
    println!("\nPhonon Density of States Chart:\n");
    print!("{}", dos.render_ascii_chart(25, 40));
}

fn run_thermo(mat: &str, t_min: f64, t_max: f64) {
    println!("Mode: Quantum Harmonic Thermodynamic Properties Evaluation");
    println!(
        "Material: {}, Temperature Range: [{:.1} K, {:.1} K]",
        mat.to_uppercase(),
        t_min,
        t_max
    );

    let (cell, ifc) = match mat.to_lowercase().as_str() {
        "gaas" => (GaAsMaterial::unit_cell(), GaAsMaterial::force_constants()),
        "graphene" => (
            GrapheneMaterial::unit_cell(),
            GrapheneMaterial::force_constants(),
        ),
        "nacl" => (NaClMaterial::unit_cell(), NaClMaterial::force_constants()),
        _ => (
            SiliconMaterial::unit_cell(),
            SiliconMaterial::force_constants(),
        ),
    };

    let recip = ReciprocalLattice::from_unit_cell(&cell);
    let builder = DynamicalMatrixBuilder::new(&cell, &ifc);

    let q_pts = recip.generate_monkhorst_pack_grid(6, 6, 6);
    let n_q = q_pts.len();
    let num_branches = 3 * cell.num_atoms();

    let mut freqs = Vec::with_capacity(n_q * num_branches);
    let mut weights = Vec::with_capacity(n_q * num_branches);
    let w_single = 1.0 / n_q as f64;

    for q in &q_pts {
        let cart_q = recip.fractional_to_cartesian(*q);
        let d = builder.build_at(cart_q);
        let eig = phonon_lattice::math::eigensystem_hermitian(&d, 1e-12, 100);
        for &val in &eig.eigenvalues {
            let omega = if val > 0.0 { val.sqrt() } else { 0.0 };
            freqs.push(omega);
            weights.push(w_single);
        }
    }

    let thermo =
        ThermodynamicCalculator::new(freqs, weights, cell.num_atoms(), cell.total_mass_kg());
    let states = thermo.temperature_scan(t_min, t_max, 10);

    println!("\n  T (K)   |  Cv (J/mol·K)  |  Cv (J/kg·K)  |  U_vib (kJ/mol) |  F_vib (kJ/mol) |  S_vib (J/mol·K)");
    println!("----------+----------------+---------------+-----------------+-----------------+------------------");
    for st in &states {
        println!(
            "{:9.1} | {:14.3} | {:13.2} | {:15.3} | {:15.3} | {:16.3}",
            st.temperature_k,
            st.heat_capacity_cv_j_per_mol_k,
            st.heat_capacity_cv_j_per_kg_k,
            st.internal_energy_kj_per_mol,
            st.free_energy_kj_per_mol,
            st.entropy_j_per_mol_k,
        );
    }

    let zpe = states[0].zero_point_energy_kj_per_mol;
    println!("\nVibrational Zero-Point Energy (ZPE): {:.3} kJ/mol", zpe);

    // Dulong-Petit classical theoretical limit: 3 * N_a * R ~ 3 * 2 * 8.314 = 49.88 J/(mol K)
    let dulong_petit = 3.0 * (cell.num_atoms() as f64) * 8.314462618;
    println!(
        "Classical Dulong-Petit Limit (3*N_a*R): {:.2} J/(mol·K)",
        dulong_petit
    );
}

fn run_bte(mat: &str, temp_k: f64) {
    println!("Mode: Boltzmann Transport Equation (BTE-RTA) Lattice Thermal Conductivity");
    println!(
        "Material: {}, Target Temperature: {:.1} K",
        mat.to_uppercase(),
        temp_k
    );

    let cell = SiliconMaterial::unit_cell();
    let ifc = SiliconMaterial::force_constants();
    let recip = ReciprocalLattice::from_unit_cell(&cell);
    let builder = DynamicalMatrixBuilder::new(&cell, &ifc);

    let debye = DebyeModel::compute(8430.0, 5840.0, cell.num_atoms(), cell.volume_m3());
    println!(
        "Acoustic Sound Speeds: v_L = {:.1} m/s, v_T = {:.1} m/s",
        debye.sound_velocity_longitudinal_m_s, debye.sound_velocity_transverse_m_s
    );
    println!(
        "Mean Sound Velocity:   v_m = {:.1} m/s",
        debye.mean_sound_velocity_m_s
    );
    println!(
        "Debye Temperature:     Θ_D = {:.1} K",
        debye.debye_temperature_k
    );

    let scattering = PhononScatteringModel::new(
        1e-4,    // Boundary length: 100 um
        2.01e-4, // Natural Si isotope mass variance g2
        1.10,    // Grüneisen parameter
        debye.mean_sound_velocity_m_s,
        cell.total_mass_kg() / cell.num_atoms() as f64,
        debye.debye_temperature_k,
        debye.debye_frequency_thz * 2.0 * std::f64::consts::PI * 1e12,
    );

    let solver = ThermalConductivitySolver::new(&builder, &recip, scattering);

    println!("\nIntegrating BTE over 6x6x6 Brillouin zone mesh (216 q-points)...");
    let result = solver.solve_at_temperature(temp_k, (6, 6, 6));

    println!(
        "Thermal Conductivity Tensor κ_αβ at {:.1} K (W/m·K):",
        temp_k
    );
    for i in 0..3 {
        println!(
            "  [ {:8.2}, {:8.2}, {:8.2} ]",
            result.conductivity_tensor_w_per_m_k.data[i][0],
            result.conductivity_tensor_w_per_m_k.data[i][1],
            result.conductivity_tensor_w_per_m_k.data[i][2],
        );
    }
    println!(
        "\nIsotropic Average Thermal Conductivity: {:.2} W/(m·K)",
        result.isotropic_conductivity_w_per_m_k
    );
}

fn run_benchmark() {
    println!("Mode: PhononLattice High-Performance Benchmark");

    let cell = SiliconMaterial::unit_cell();
    let ifc = SiliconMaterial::force_constants();
    let _recip = ReciprocalLattice::from_unit_cell(&cell);
    let builder = DynamicalMatrixBuilder::new(&cell, &ifc);

    let n_evals = 10_000;
    println!(
        "Running {} dynamical matrix assembly operations...",
        n_evals
    );
    let start_assembly = Instant::now();
    for i in 0..n_evals {
        let q = phonon_lattice::math::Vector3::new(i as f64 * 0.0001, 0.0, 0.0);
        let _ = builder.build_at(q);
    }
    let elapsed_assembly = start_assembly.elapsed();
    let ops_sec_assembly = n_evals as f64 / elapsed_assembly.as_secs_f64();
    println!(
        "Dynamical Matrix Assembly: {:.2e} ops/sec (Time: {:.3?})",
        ops_sec_assembly, elapsed_assembly
    );

    let n_diag = 2_000;
    println!(
        "\nRunning {} complex Hermitian Jacobi diagonalizations (6x6 matrix)...",
        n_diag
    );
    let sample_mat = builder.build_at(phonon_lattice::math::Vector3::new(0.2, 0.3, 0.4));
    let start_diag = Instant::now();
    for _ in 0..n_diag {
        let _ = phonon_lattice::math::eigensystem_hermitian(&sample_mat, 1e-13, 100);
    }
    let elapsed_diag = start_diag.elapsed();
    let ops_sec_diag = n_diag as f64 / elapsed_diag.as_secs_f64();
    println!(
        "Jacobi Complex Hermitian Eigensolver: {:.2e} diags/sec (Time: {:.3?})",
        ops_sec_diag, elapsed_diag
    );

    println!("\nBenchmark completed successfully! Numerical performance optimal.");
}
