use crate::thermodynamics::harmonic::BOLTZMANN_CONSTANT;
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ScatteringRates {
    pub boundary_rate_s_inv: f64,
    pub isotope_rate_s_inv: f64,
    pub umklapp_rate_s_inv: f64,
    pub total_rate_s_inv: f64,
    pub relaxation_time_s: f64,
    pub mean_free_path_nm: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PhononScatteringModel {
    /// Effective Casimir grain boundary length in meters (e.g. 1e-4 m = 100 um)
    pub boundary_length_m: f64,
    /// Isotope mass variance factor g2 (e.g. 2.01e-4 for natural Silicon)
    pub isotope_variance_g2: f64,
    /// Mode Gruneisen parameter gamma (typically ~1.0 to 1.5)
    pub gruneisen_parameter: f64,
    /// Average acoustic sound speed in m/s
    pub sound_speed_m_s: f64,
    /// Average atomic mass in kg
    pub average_mass_kg: f64,
    /// Debye temperature in K
    pub debye_temperature_k: f64,
    /// Debye frequency in rad/s
    pub debye_frequency_rad_s: f64,
}

impl PhononScatteringModel {
    pub fn new(
        boundary_length_m: f64,
        isotope_variance_g2: f64,
        gruneisen_parameter: f64,
        sound_speed_m_s: f64,
        average_mass_kg: f64,
        debye_temperature_k: f64,
        debye_frequency_rad_s: f64,
    ) -> Self {
        Self {
            boundary_length_m,
            isotope_variance_g2,
            gruneisen_parameter,
            sound_speed_m_s,
            average_mass_kg,
            debye_temperature_k,
            debye_frequency_rad_s,
        }
    }

    /// Evaluates scattering rates and relaxation time tau(omega, T) for a mode with frequency omega (rad/s)
    /// and group velocity v_g (m/s) at temperature T (K).
    pub fn compute_scattering_rates(
        &self,
        omega_rad_s: f64,
        v_group_m_s: f64,
        temp_k: f64,
    ) -> ScatteringRates {
        if omega_rad_s <= 1e6 || v_group_m_s <= 1.0 {
            return ScatteringRates {
                boundary_rate_s_inv: 0.0,
                isotope_rate_s_inv: 0.0,
                umklapp_rate_s_inv: 0.0,
                total_rate_s_inv: 1e15, // Extremely high damping / near zero lifetime
                relaxation_time_s: 0.0,
                mean_free_path_nm: 0.0,
            };
        }

        let temp = temp_k.max(0.1);

        // 1. Boundary scattering rate: 1 / tau_B = v_g / L
        let rate_boundary = v_group_m_s / self.boundary_length_m.max(1e-9);

        // 2. Isotope scattering rate (Klemens-Tamura formula):
        // 1 / tau_iso = (V_0 / (4 * pi * v_s^3)) * g2 * omega^4
        // Approximated using Debye density of states: 1 / tau_iso = (pi / 2) * g2 * (omega^4 / omega_D^3)
        let rate_iso = (std::f64::consts::PI * 0.5)
            * self.isotope_variance_g2
            * (omega_rad_s.powi(4) / self.debye_frequency_rad_s.powi(3).max(1e10));

        // 3. Three-phonon Umklapp scattering rate (Klemens-Callaway / Slack):
        // 1 / tau_U = 2 * gamma^2 * (k_B * T / (M * v_s^2)) * (omega^2 / omega_D) * exp(-Theta_D / (3 * T))
        let mass_v2 = self.average_mass_kg * self.sound_speed_m_s.powi(2);
        let thermal_factor = (BOLTZMANN_CONSTANT * temp) / mass_v2.max(1e-30);
        let exp_umklapp = (-self.debye_temperature_k / (3.0 * temp)).exp();

        let rate_umklapp = 2.0
            * self.gruneisen_parameter.powi(2)
            * thermal_factor
            * (omega_rad_s.powi(2) / self.debye_frequency_rad_s.max(1e6))
            * exp_umklapp;

        // Total scattering rate (Matthiessen's rule)
        let rate_total = rate_boundary + rate_iso + rate_umklapp;
        let tau = 1.0 / rate_total;
        let mfp_nm = (v_group_m_s * tau) * 1e9; // in nm

        ScatteringRates {
            boundary_rate_s_inv: rate_boundary,
            isotope_rate_s_inv: rate_iso,
            umklapp_rate_s_inv: rate_umklapp,
            total_rate_s_inv: rate_total,
            relaxation_time_s: tau,
            mean_free_path_nm: mfp_nm,
        }
    }
}
