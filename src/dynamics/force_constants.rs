use crate::math::Matrix3;
use serde::{Deserialize, Serialize};

pub const EV_PER_ANGSTROM2_TO_N_PER_M: f64 = 16.02176634;

/// 3x3 force constant matrix Phi_alpha_beta(0j, Rj')
pub type ForceConstantTensor = Matrix3;

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ForceConstantEntry {
    /// Index of atom in reference cell (0, 0, 0)
    pub atom_idx_0: usize,
    /// Lattice image indices (n1, n2, n3)
    pub cell_image: (i32, i32, i32),
    /// Index of atom in target cell R
    pub atom_idx_r: usize,
    /// 3x3 tensor in N/m (J/m^2)
    pub tensor_n_per_m: Matrix3,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct InteratomicForceConstants {
    pub num_atoms: usize,
    pub entries: Vec<ForceConstantEntry>,
}

impl InteratomicForceConstants {
    pub fn new(num_atoms: usize) -> Self {
        Self {
            num_atoms,
            entries: Vec::new(),
        }
    }

    /// Adds a force constant tensor between atom j in cell (0,0,0) and atom j' in cell (n1, n2, n3).
    /// Input tensor can be provided in eV/Angstrom^2 (converted automatically) or directly in N/m.
    pub fn add_entry_ev_per_angstrom2(
        &mut self,
        atom_idx_0: usize,
        cell_image: (i32, i32, i32),
        atom_idx_r: usize,
        tensor_ev_ang2: Matrix3,
    ) {
        let tensor_n_per_m = tensor_ev_ang2.mul_scalar(EV_PER_ANGSTROM2_TO_N_PER_M);
        self.entries.push(ForceConstantEntry {
            atom_idx_0,
            cell_image,
            atom_idx_r,
            tensor_n_per_m,
        });
    }

    pub fn add_entry_n_per_m(
        &mut self,
        atom_idx_0: usize,
        cell_image: (i32, i32, i32),
        atom_idx_r: usize,
        tensor_n_per_m: Matrix3,
    ) {
        self.entries.push(ForceConstantEntry {
            atom_idx_0,
            cell_image,
            atom_idx_r,
            tensor_n_per_m,
        });
    }

    /// Enforces the Acoustic Sum Rule (ASR):
    /// Phi(0j, 0j) = - sum_{(R, j') != (0, j)} Phi(0j, Rj')
    /// Ensures that rigid translations incur zero restoring force and acoustic phonon modes at Gamma strictly equal 0.
    pub fn enforce_acoustic_sum_rule(&mut self) {
        for j in 0..self.num_atoms {
            let mut sum_other = Matrix3::ZERO;

            for entry in &self.entries {
                if entry.atom_idx_0 == j {
                    let is_self = entry.atom_idx_r == j && entry.cell_image == (0, 0, 0);
                    if !is_self {
                        sum_other = sum_other + entry.tensor_n_per_m;
                    }
                }
            }

            let required_self = Matrix3::ZERO - sum_other;

            // Find existing (0j, 0j) entry or append it
            let mut found = false;
            for entry in &mut self.entries {
                if entry.atom_idx_0 == j && entry.atom_idx_r == j && entry.cell_image == (0, 0, 0) {
                    entry.tensor_n_per_m = required_self;
                    found = true;
                    break;
                }
            }

            if !found {
                self.entries.push(ForceConstantEntry {
                    atom_idx_0: j,
                    cell_image: (0, 0, 0),
                    atom_idx_r: j,
                    tensor_n_per_m: required_self,
                });
            }
        }
    }
}
