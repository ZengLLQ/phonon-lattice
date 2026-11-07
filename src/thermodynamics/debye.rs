use super::harmonic::{BOLTZMANN_CONSTANT, HBAR};
use serde::{Deserialize, Serialize};
use std::f64::consts::PI;

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct DebyeProperties {
    pub sound_velocity_longitudinal_m_s: f64,
    pub sound_velocity_transverse_m_s: f64,
    pub mean_sound_velocity_m_s: f64,
    pub debye_frequency_thz: f64,
    pub debye_temperature_k: f64,
}

pub struct DebyeModel;

impl DebyeModel {
    /// Computes the isotropic mean sound velocity v_m:
    /// 3 / v_m^3 = 1 / v_L^3 + 2 / v_T^3
    pub fn mean_sound_velocity(v_l: f64, v_t: f64) -> f64 {
        let inv_v3 = (1.0 / v_l.powi(3) + 2.0 / v_t.powi(3)) / 3.0;
        (1.0 / inv_v3).cbrt()
    }

    /// Evaluates Debye properties from acoustic sound speeds and unit cell volume / atom count
    pub fn compute(
        v_l: f64,
        v_t: f64,
        num_atoms_per_cell: usize,
        cell_volume_m3: f64,
    ) -> DebyeProperties {
        let v_m = Self::mean_sound_velocity(v_l, v_t);
        let atom_density = num_atoms_per_cell as f64 / cell_volume_m3;

        // Debye cut-off wavevector: q_D = (6 * pi^2 * n)^(1/3)
        let q_debye = (6.0 * PI * PI * atom_density).cbrt();

        // Debye cutoff angular frequency: omega_D = v_m * q_D
        let omega_debye = v_m * q_debye;
        let debye_freq_thz = omega_debye / (2.0 * PI * 1e12);

        // Debye temperature: Theta_D = hbar * omega_D / k_B
        let theta_debye = (HBAR * omega_debye) / BOLTZMANN_CONSTANT;

        DebyeProperties {
            sound_velocity_longitudinal_m_s: v_l,
            sound_velocity_transverse_m_s: v_t,
            mean_sound_velocity_m_s: v_m,
            debye_frequency_thz: debye_freq_thz,
            debye_temperature_k: theta_debye,
        }
    }

    /// Low temperature Debye T^3 heat capacity approximation in J / (mol K):
    /// C_v(T) = (12 * pi^4 / 5) * N_a * k_B * (T / Theta_D)^3
    pub fn low_temp_cv(temp_k: f64, theta_d: f64, num_atoms_per_cell: usize) -> f64 {
        let r = temp_k / theta_d;
        if r > 0.15 {
            // Beyond valid T^3 range
            return 0.0;
        }
        let prefactor = (12.0 * PI.powi(4) / 5.0)
            * (num_atoms_per_cell as f64)
            * (super::harmonic::AVOGADRO_CONSTANT * BOLTZMANN_CONSTANT);
        prefactor * r.powi(3)
    }
}
