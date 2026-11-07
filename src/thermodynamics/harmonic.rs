use serde::{Deserialize, Serialize};

pub const BOLTZMANN_CONSTANT: f64 = 1.380649e-23; // J/K
pub const HBAR: f64 = 1.054571817e-34; // J s
pub const AVOGADRO_CONSTANT: f64 = 6.02214076e23; // mol^-1
pub const JOULE_TO_EV: f64 = 6.241509074e18; // eV/J

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ThermodynamicState {
    pub temperature_k: f64,
    /// Constant-volume heat capacity in J / (mol K)
    pub heat_capacity_cv_j_per_mol_k: f64,
    /// Constant-volume heat capacity in J / (kg K)
    pub heat_capacity_cv_j_per_kg_k: f64,
    /// Vibrational internal energy in kJ / mol
    pub internal_energy_kj_per_mol: f64,
    /// Vibrational Helmholtz free energy in kJ / mol
    pub free_energy_kj_per_mol: f64,
    /// Vibrational entropy in J / (mol K)
    pub entropy_j_per_mol_k: f64,
    /// Zero-point vibrational energy in kJ / mol
    pub zero_point_energy_kj_per_mol: f64,
}

pub struct ThermodynamicCalculator {
    /// Distinct phonon frequencies collected over reciprocal mesh in rad/s
    pub frequencies_rad_s: Vec<f64>,
    /// Statistical weight of each frequency (sum of weights = 1.0)
    pub weights: Vec<f64>,
    /// Number of atoms per primitive cell
    pub num_atoms: usize,
    /// Total mass of unit cell in kg
    pub cell_mass_kg: f64,
}

impl ThermodynamicCalculator {
    pub fn new(
        frequencies_rad_s: Vec<f64>,
        weights: Vec<f64>,
        num_atoms: usize,
        cell_mass_kg: f64,
    ) -> Self {
        assert_eq!(
            frequencies_rad_s.len(),
            weights.len(),
            "Frequency and weight counts must match"
        );
        Self {
            frequencies_rad_s,
            weights,
            num_atoms,
            cell_mass_kg,
        }
    }

    /// Evaluates thermodynamic state at temperature T (Kelvin)
    pub fn evaluate_at_temperature(&self, temp_k: f64) -> ThermodynamicState {
        let temp = temp_k.max(1e-6); // Prevent division by zero
        let kt = BOLTZMANN_CONSTANT * temp;

        let mut zpe_cell = 0.0;
        let mut u_cell = 0.0;
        let mut f_cell = 0.0;
        let mut s_cell = 0.0;
        let mut cv_cell = 0.0;

        for (i, &omega) in self.frequencies_rad_s.iter().enumerate() {
            if omega <= 1e6 {
                // Skip zero acoustic frequency at Gamma
                continue;
            }

            let w = self.weights[i];
            let hw = HBAR * omega;
            let x = hw / kt;

            // Zero-point energy per mode
            zpe_cell += w * 0.5 * hw;

            if x > 100.0 {
                // High energy / low T limit: zero occupation
                u_cell += w * 0.5 * hw;
                f_cell += w * 0.5 * hw;
                // cv_cell and s_cell -> 0
            } else if x < 1e-4 {
                // Classical equipartition limit: T -> infinity
                let exp_x = x.exp();
                let n_occ = 1.0 / (exp_x - 1.0);
                u_cell += w * (0.5 + n_occ) * hw;
                f_cell += w * (0.5 * hw + kt * (1.0 - (-x).exp()).ln());
                s_cell += w * BOLTZMANN_CONSTANT * (n_occ * x - (1.0 - (-x).exp()).ln());
                // Dulong-Petit contribution: k_B * (x^2 * exp(x) / (exp(x)-1)^2) ~ k_B
                cv_cell += w * BOLTZMANN_CONSTANT;
            } else {
                let exp_x = x.exp();
                let n_occ = 1.0 / (exp_x - 1.0);

                u_cell += w * (0.5 + n_occ) * hw;
                f_cell += w * (0.5 * hw + kt * (1.0 - (-x).exp()).ln());
                s_cell += w * BOLTZMANN_CONSTANT * (n_occ * x - (1.0 - (-x).exp()).ln());

                let cv_contrib = BOLTZMANN_CONSTANT * (x * x * exp_x) / ((exp_x - 1.0).powi(2));
                cv_cell += w * cv_contrib;
            }
        }

        // Convert per-cell values to molar quantities (per mole of formula units)
        let zpe_kj_mol = zpe_cell * AVOGADRO_CONSTANT * 1e-3;
        let u_kj_mol = u_cell * AVOGADRO_CONSTANT * 1e-3;
        let f_kj_mol = f_cell * AVOGADRO_CONSTANT * 1e-3;
        let s_j_mol_k = s_cell * AVOGADRO_CONSTANT;
        let cv_j_mol_k = cv_cell * AVOGADRO_CONSTANT;

        // Specific heat capacity per kg
        let mass_per_mol_kg = self.cell_mass_kg * AVOGADRO_CONSTANT;
        let cv_j_kg_k = cv_j_mol_k / mass_per_mol_kg;

        ThermodynamicState {
            temperature_k: temp_k,
            heat_capacity_cv_j_per_mol_k: cv_j_mol_k,
            heat_capacity_cv_j_per_kg_k: cv_j_kg_k,
            internal_energy_kj_per_mol: u_kj_mol,
            free_energy_kj_per_mol: f_kj_mol,
            entropy_j_per_mol_k: s_j_mol_k,
            zero_point_energy_kj_per_mol: zpe_kj_mol,
        }
    }

    /// Evaluates thermodynamic curve over temperature range [t_min, t_max] with n_points
    pub fn temperature_scan(
        &self,
        t_min: f64,
        t_max: f64,
        n_points: usize,
    ) -> Vec<ThermodynamicState> {
        assert!(n_points >= 2, "Must specify at least 2 temperature points");
        let dt = (t_max - t_min) / (n_points as f64 - 1.0);
        (0..n_points)
            .map(|i| {
                let t = t_min + i as f64 * dt;
                self.evaluate_at_temperature(t)
            })
            .collect()
    }
}
