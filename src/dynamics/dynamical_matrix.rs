use super::force_constants::InteratomicForceConstants;
use crate::crystal::UnitCell;
use crate::math::Vector3;
use num_complex::Complex64;

pub const VACUUM_PERMITTIVITY: f64 = 8.8541878128e-12; // F/m = C^2 / (N m^2)
pub const ELEMENTARY_CHARGE: f64 = 1.602176634e-19; // C

pub struct DynamicalMatrixBuilder<'a> {
    pub cell: &'a UnitCell,
    pub ifc: &'a InteratomicForceConstants,
}

impl<'a> DynamicalMatrixBuilder<'a> {
    pub fn new(cell: &'a UnitCell, ifc: &'a InteratomicForceConstants) -> Self {
        assert_eq!(
            cell.num_atoms(),
            ifc.num_atoms,
            "Cell atom count must match force constant atom count"
        );
        Self { cell, ifc }
    }

    /// Evaluates the (3 * N_atoms) x (3 * N_atoms) complex Hermitian dynamical matrix D(q) in s^-2.
    /// q_cart is in Angstroms^-1.
    pub fn build_at(&self, q_cart: Vector3) -> Vec<Vec<Complex64>> {
        let n_atoms = self.cell.num_atoms();
        let dim = 3 * n_atoms;
        let mut d = vec![vec![Complex64::new(0.0, 0.0); dim]; dim];

        // 1. Harmonic force constant contribution:
        // D_ab(j, j'; q) = 1/sqrt(M_j * M_j') * sum_R Phi_ab(0j, Rj') * exp(i * q . R)
        for entry in &self.ifc.entries {
            let j = entry.atom_idx_0;
            let j_prime = entry.atom_idx_r;
            let (n1, n2, n3) = entry.cell_image;

            let m_j = self.cell.atoms[j].mass_kg();
            let m_j_prime = self.cell.atoms[j_prime].mass_kg();
            let mass_factor = 1.0 / (m_j * m_j_prime).sqrt();

            // Lattice vector R in Angstroms
            let r_vec = self.cell.a1 * (n1 as f64)
                + self.cell.a2 * (n2 as f64)
                + self.cell.a3 * (n3 as f64);

            // Phase factor q . R (dimensionless because q is in Angstrom^-1 and R is in Angstrom)
            let q_dot_r = q_cart.dot(&r_vec);
            let phase = Complex64::from_polar(1.0, q_dot_r);

            let tensor = entry.tensor_n_per_m; // in N/m

            for alpha in 0..3 {
                for beta in 0..3 {
                    let row = 3 * j + alpha;
                    let col = 3 * j_prime + beta;
                    let phi_val = tensor.data[alpha][beta] * mass_factor; // in s^-2
                    d[row][col] += phase * phi_val;
                }
            }
        }

        // 2. Non-analytical correction for polar crystals (LO-TO splitting) as q -> 0
        let q_norm = q_cart.norm();
        if q_norm > 1e-6 {
            let q_dir = q_cart.normalize();
            let vol_m3 = self.cell.volume_m3();
            let eps = self.cell.dielectric_tensor;

            // Denominator: q_dir . epsilon . q_dir
            let eps_q = eps.mul_vec(q_dir);
            let denom = q_dir.dot(&eps_q);

            if denom > 1e-6 {
                // Prefactor: 1 / (epsilon_0 * Omega_0)
                let prefactor = 1.0 / (VACUUM_PERMITTIVITY * vol_m3 * denom);

                // D_NA(j, j') = prefactor / sqrt(M_j M_j') * (q_dir . Z_j*)_alpha * (q_dir . Z_j'* )_beta
                for j in 0..n_atoms {
                    if let Some(z_j) = self.cell.atoms[j].born_charge {
                        // q_dir . Z_j (vector product in units of elementary charge e)
                        let q_z_j = z_j.transpose().mul_vec(q_dir) * ELEMENTARY_CHARGE;
                        let m_j = self.cell.atoms[j].mass_kg();

                        for j_prime in 0..n_atoms {
                            if let Some(z_j_prime) = self.cell.atoms[j_prime].born_charge {
                                let q_z_jp =
                                    z_j_prime.transpose().mul_vec(q_dir) * ELEMENTARY_CHARGE;
                                let m_jp = self.cell.atoms[j_prime].mass_kg();
                                let m_factor = 1.0 / (m_j * m_jp).sqrt();

                                for alpha in 0..3 {
                                    for beta in 0..3 {
                                        let q_comp_j = match alpha {
                                            0 => q_z_j.x,
                                            1 => q_z_j.y,
                                            _ => q_z_j.z,
                                        };
                                        let q_comp_jp = match beta {
                                            0 => q_z_jp.x,
                                            1 => q_z_jp.y,
                                            _ => q_z_jp.z,
                                        };

                                        let val_na = prefactor * m_factor * q_comp_j * q_comp_jp;
                                        let row = 3 * j + alpha;
                                        let col = 3 * j_prime + beta;
                                        d[row][col] += Complex64::new(val_na, 0.0);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        // 3. Enforce Hermiticity: D_dagger = D
        for i in 0..dim {
            d[i][i] = Complex64::new(d[i][i].re, 0.0);
            for j in (i + 1)..dim {
                let avg = (d[i][j] + d[j][i].conj()) * 0.5;
                d[i][j] = avg;
                d[j][i] = avg.conj();
            }
        }

        d
    }

    /// Computes analytical group velocities v_g,s(q) = grad_q omega_s(q) using central finite differences:
    /// dD/dq_alpha ~ (D(q + dq_alpha) - D(q - dq_alpha)) / (2 * dq)
    /// v_g,s = 1 / (2 * omega_s) * <e_s | dD/dq | e_s>
    pub fn compute_group_velocities(
        &self,
        q_cart: Vector3,
        frequencies_rad_s: &[f64],
        eigenvectors: &[Vec<Complex64>],
    ) -> Vec<Vector3> {
        let dq = 1e-4; // in Angstroms^-1 (10^6 m^-1)
        let dq_m = dq * 1e10; // in m^-1

        let num_modes = frequencies_rad_s.len();
        let mut group_velocities = Vec::with_capacity(num_modes);

        let d_plus_x = self.build_at(q_cart + Vector3::new(dq, 0.0, 0.0));
        let d_minus_x = self.build_at(q_cart - Vector3::new(dq, 0.0, 0.0));

        let d_plus_y = self.build_at(q_cart + Vector3::new(0.0, dq, 0.0));
        let d_minus_y = self.build_at(q_cart - Vector3::new(0.0, dq, 0.0));

        let d_plus_z = self.build_at(q_cart + Vector3::new(0.0, 0.0, dq));
        let d_minus_z = self.build_at(q_cart - Vector3::new(0.0, 0.0, dq));

        let dim = 3 * self.cell.num_atoms();

        for s in 0..num_modes {
            let omega = frequencies_rad_s[s];
            if omega < 1e6 {
                // Near Gamma acoustic mode
                group_velocities.push(Vector3::ZERO);
                continue;
            }

            let e_s = &eigenvectors[s];

            let mut grad_d_x = Complex64::new(0.0, 0.0);
            let mut grad_d_y = Complex64::new(0.0, 0.0);
            let mut grad_d_z = Complex64::new(0.0, 0.0);

            for i in 0..dim {
                for j in 0..dim {
                    let factor = e_s[i].conj() * e_s[j];
                    let diff_x = (d_plus_x[i][j] - d_minus_x[i][j]) / (2.0 * dq_m);
                    let diff_y = (d_plus_y[i][j] - d_minus_y[i][j]) / (2.0 * dq_m);
                    let diff_z = (d_plus_z[i][j] - d_minus_z[i][j]) / (2.0 * dq_m);

                    grad_d_x += factor * diff_x;
                    grad_d_y += factor * diff_y;
                    grad_d_z += factor * diff_z;
                }
            }

            let vx = grad_d_x.re / (2.0 * omega);
            let vy = grad_d_y.re / (2.0 * omega);
            let vz = grad_d_z.re / (2.0 * omega);

            group_velocities.push(Vector3::new(vx, vy, vz));
        }

        group_velocities
    }
}
