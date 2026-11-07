use super::scattering::PhononScatteringModel;
use crate::crystal::ReciprocalLattice;
use crate::dynamics::DynamicalMatrixBuilder;
use crate::math::{eigensystem_hermitian, Matrix3};
use crate::thermodynamics::harmonic::{BOLTZMANN_CONSTANT, HBAR};
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ThermalConductivityResult {
    pub temperature_k: f64,
    /// Physical symmetrized 3x3 thermal conductivity tensor in W / (m K)
    pub conductivity_tensor_w_per_m_k: Matrix3,
    /// Raw unsymmetrized tensor evaluated directly from reciprocal grid sum
    pub raw_tensor_w_per_m_k: Matrix3,
    /// Isotropic average conductivity: tr(kappa) / 3
    pub isotropic_conductivity_w_per_m_k: f64,
    /// Total number of phonon modes integrated
    pub total_modes: usize,
}

pub struct ThermalConductivitySolver<'a> {
    pub builder: &'a DynamicalMatrixBuilder<'a>,
    pub reciprocal: &'a ReciprocalLattice,
    pub scattering: PhononScatteringModel,
}

impl<'a> ThermalConductivitySolver<'a> {
    pub fn new(
        builder: &'a DynamicalMatrixBuilder<'a>,
        reciprocal: &'a ReciprocalLattice,
        scattering: PhononScatteringModel,
    ) -> Self {
        Self {
            builder,
            reciprocal,
            scattering,
        }
    }

    /// Solves the Boltzmann Transport Equation for lattice thermal conductivity at temperature T (K)
    /// across an (n1 x n2 x n3) Monkhorst-Pack Brillouin zone grid.
    pub fn solve_at_temperature(
        &self,
        temp_k: f64,
        grid_dims: (usize, usize, usize),
    ) -> ThermalConductivityResult {
        let temp = temp_k.max(0.1);
        let kt = BOLTZMANN_CONSTANT * temp;
        let cell_vol_m3 = self.builder.cell.volume_m3();

        let q_points =
            self.reciprocal
                .generate_monkhorst_pack_grid(grid_dims.0, grid_dims.1, grid_dims.2);
        let n_q = q_points.len();
        let n_atoms = self.builder.cell.num_atoms();
        let num_branches = 3 * n_atoms;

        let mut kappa_tensor = [[0.0; 3]; 3];
        let mut valid_modes = 0;

        for frac_q in &q_points {
            let cart_q = self.reciprocal.fractional_to_cartesian(*frac_q);
            let d_mat = self.builder.build_at(cart_q);
            let eigen = eigensystem_hermitian(&d_mat, 1e-12, 100);

            let mut freqs_rad_s = Vec::with_capacity(num_branches);
            for &val in &eigen.eigenvalues {
                freqs_rad_s.push(if val > 0.0 { val.sqrt() } else { 0.0 });
            }

            let group_vels =
                self.builder
                    .compute_group_velocities(cart_q, &freqs_rad_s, &eigen.eigenvectors);

            for s in 0..num_branches {
                let omega = freqs_rad_s[s];
                if omega <= 1e6 {
                    continue;
                }

                let v_vec = group_vels[s];
                let v_norm = v_vec.norm();
                if v_norm < 1.0 {
                    continue;
                }

                // Mode heat capacity: C_v(q, s)
                let hw = HBAR * omega;
                let x = hw / kt;
                let c_mode = if x > 80.0 {
                    0.0
                } else if x < 1e-4 {
                    BOLTZMANN_CONSTANT
                } else {
                    let exp_x = x.exp();
                    BOLTZMANN_CONSTANT * (x * x * exp_x) / ((exp_x - 1.0).powi(2))
                };

                if c_mode < 1e-35 {
                    continue;
                }

                // Mode relaxation time tau(q, s)
                let rates = self
                    .scattering
                    .compute_scattering_rates(omega, v_norm, temp);
                let tau = rates.relaxation_time_s;

                // kappa_ab += C_mode * v_a * v_b * tau
                let v_comps = [v_vec.x, v_vec.y, v_vec.z];
                for a in 0..3 {
                    for b in 0..3 {
                        kappa_tensor[a][b] += c_mode * v_comps[a] * v_comps[b] * tau;
                    }
                }

                valid_modes += 1;
            }
        }

        // Normalization: 1 / (V_cell * N_q)
        let norm_factor = 1.0 / (cell_vol_m3 * n_q as f64);
        for a in 0..3 {
            for b in 0..3 {
                kappa_tensor[a][b] *= norm_factor;
            }
        }

        let raw_matrix = Matrix3::new(kappa_tensor);
        let isotropic = raw_matrix.trace() / 3.0;

        // Apply cubic point group symmetry projection: kappa_symm = tr(kappa)/3 * I
        let symmetrized = Matrix3::IDENTITY.mul_scalar(isotropic);

        ThermalConductivityResult {
            temperature_k: temp,
            conductivity_tensor_w_per_m_k: symmetrized,
            raw_tensor_w_per_m_k: raw_matrix,
            isotropic_conductivity_w_per_m_k: isotropic,
            total_modes: valid_modes,
        }
    }

    /// Evaluates thermal conductivity temperature curve from t_min to t_max
    pub fn temperature_scan(
        &self,
        t_min: f64,
        t_max: f64,
        n_points: usize,
        grid_dims: (usize, usize, usize),
    ) -> Vec<ThermalConductivityResult> {
        let dt = (t_max - t_min) / (n_points as f64 - 1.0);
        (0..n_points)
            .map(|i| {
                let t = t_min + i as f64 * dt;
                self.solve_at_temperature(t, grid_dims)
            })
            .collect()
    }
}
