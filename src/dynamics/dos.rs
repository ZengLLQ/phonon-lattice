use super::dispersion::RAD_S_TO_THZ;
use super::dynamical_matrix::DynamicalMatrixBuilder;
use crate::crystal::ReciprocalLattice;
use crate::math::{eigensystem_hermitian, gaussian_delta, integrate_trapezoidal};
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DensityOfStates {
    /// Frequency mesh in THz
    pub frequencies_thz: Vec<f64>,
    /// Total phonon density of states in states/THz
    pub total_dos: Vec<f64>,
    /// Projected density of states for each atom index [atom_idx][freq_idx]
    pub projected_dos: Vec<Vec<f64>>,
    /// Integrated DOS sum (should equal 3 * N_atoms)
    pub integrated_dos: f64,
}

impl DensityOfStates {
    /// Renders an ASCII horizontal bar chart preview of the DOS
    pub fn render_ascii_chart(&self, num_bins: usize, max_bar_width: usize) -> String {
        let n = self.frequencies_thz.len();
        if n == 0 {
            return String::new();
        }

        let max_val = self.total_dos.iter().cloned().fold(0.0, f64::max);
        if max_val < 1e-12 {
            return "Empty DOS".to_string();
        }

        let step = (n / num_bins).max(1);
        let mut out = String::new();
        out.push_str("Frequency (THz) | Density of States (DOS)\n");
        out.push_str("----------------+--------------------------------------------------\n");

        for i in (0..n).step_by(step) {
            let freq = self.frequencies_thz[i];
            let val = self.total_dos[i];
            let bar_len = ((val / max_val) * max_bar_width as f64).round() as usize;
            let bar: String = "█".repeat(bar_len);
            out.push_str(&format!("{:15.2} | {:5.2} {}\n", freq, val, bar));
        }

        out
    }
}

pub struct PhononDosCalculator<'a> {
    pub builder: &'a DynamicalMatrixBuilder<'a>,
    pub reciprocal: &'a ReciprocalLattice,
}

impl<'a> PhononDosCalculator<'a> {
    pub fn new(builder: &'a DynamicalMatrixBuilder<'a>, reciprocal: &'a ReciprocalLattice) -> Self {
        Self {
            builder,
            reciprocal,
        }
    }

    /// Computes Total DOS and Atom-Projected DOS over a uniform Monkhorst-Pack grid.
    /// grid_dims: (n1, n2, n3)
    /// sigma_thz: Gaussian broadening parameter in THz (typically 0.1 - 0.2 THz)
    /// num_energy_points: number of frequency bins
    pub fn compute(
        &self,
        grid_dims: (usize, usize, usize),
        sigma_thz: f64,
        num_energy_points: usize,
    ) -> DensityOfStates {
        let q_points =
            self.reciprocal
                .generate_monkhorst_pack_grid(grid_dims.0, grid_dims.1, grid_dims.2);
        let n_q = q_points.len();
        let n_atoms = self.builder.cell.num_atoms();
        let num_branches = 3 * n_atoms;

        // Collect all eigenvalues and atomic projection weights
        // all_modes: Vec<(frequency_thz, Vec<atom_projection_weight>)>
        let mut all_modes = Vec::with_capacity(n_q * num_branches);
        let mut max_freq = 0.0;

        for frac_q in &q_points {
            let cart_q = self.reciprocal.fractional_to_cartesian(*frac_q);
            let d_mat = self.builder.build_at(cart_q);
            let eigen = eigensystem_hermitian(&d_mat, 1e-12, 100);

            for (s, &val) in eigen.eigenvalues.iter().enumerate() {
                let freq = if val > 0.0 {
                    val.sqrt() * RAD_S_TO_THZ
                } else {
                    0.0
                };

                if freq > max_freq {
                    max_freq = freq;
                }

                let e_s = &eigen.eigenvectors[s];
                let mut atom_weights = Vec::with_capacity(n_atoms);
                for j in 0..n_atoms {
                    let mut proj = 0.0;
                    for alpha in 0..3 {
                        proj += e_s[3 * j + alpha].norm_sqr();
                    }
                    atom_weights.push(proj);
                }

                all_modes.push((freq, atom_weights));
            }
        }

        // Energy grid from 0.0 to 1.1 * max_freq
        let upper_freq = (max_freq * 1.1).max(1.0);
        let d_freq = upper_freq / (num_energy_points as f64 - 1.0);
        let mut freqs_thz = Vec::with_capacity(num_energy_points);
        let mut total_dos = vec![0.0; num_energy_points];
        let mut projected_dos = vec![vec![0.0; num_energy_points]; n_atoms];

        let weight = 1.0 / (n_q as f64);

        for i in 0..num_energy_points {
            let omega = i as f64 * d_freq;
            freqs_thz.push(omega);

            for (mode_freq, atom_weights) in &all_modes {
                let delta = gaussian_delta(omega - mode_freq, sigma_thz);
                let contrib = weight * delta;
                total_dos[i] += contrib;

                for j in 0..n_atoms {
                    projected_dos[j][i] += contrib * atom_weights[j];
                }
            }
        }

        let integrated = integrate_trapezoidal(&total_dos, d_freq);

        DensityOfStates {
            frequencies_thz: freqs_thz,
            total_dos,
            projected_dos,
            integrated_dos: integrated,
        }
    }
}
