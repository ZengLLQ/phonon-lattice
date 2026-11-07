use crate::crystal::{Atom, HighSymmetryPath, ReciprocalLattice, UnitCell};
use crate::dynamics::force_constants::InteratomicForceConstants;
use crate::math::{Matrix3, Vector3};

pub struct SiliconMaterial;

impl SiliconMaterial {
    /// Diamond Silicon lattice constant a_0 = 5.431 Angstroms
    pub const LATTICE_CONSTANT: f64 = 5.431;
    pub const MASS_AMU: f64 = 28.0855;

    pub fn unit_cell() -> UnitCell {
        let a = Self::LATTICE_CONSTANT * 0.5;
        // FCC primitive vectors
        let a1 = Vector3::new(0.0, a, a);
        let a2 = Vector3::new(a, 0.0, a);
        let a3 = Vector3::new(a, a, 0.0);

        let atoms = vec![
            Atom::new("Si", Self::MASS_AMU, Vector3::new(0.0, 0.0, 0.0)),
            Atom::new("Si", Self::MASS_AMU, Vector3::new(0.25, 0.25, 0.25)),
        ];

        // Silicon is non-polar, epsilon_inf ~ 11.9
        UnitCell::new(a1, a2, a3, atoms).with_dielectric_tensor(Matrix3::IDENTITY.mul_scalar(11.9))
    }

    pub fn force_constants() -> InteratomicForceConstants {
        let mut ifc = InteratomicForceConstants::new(2);

        // Nearest neighbor bond-stretching and bond-bending parameters (in eV/Angstrom^2)
        // Yields experimental Silicon optical phonon frequency at Gamma: 15.6 THz (520 cm^-1)
        let alpha = 6.28;
        let beta = 2.10;

        // 4 nearest neighbors of atom 0 (atom 1 in various cells)
        let nn_bonds = [
            ((0, 0, 0), Vector3::new(1.0, 1.0, 1.0).normalize()),
            ((-1, 0, 0), Vector3::new(-1.0, -1.0, 1.0).normalize()),
            ((0, -1, 0), Vector3::new(1.0, -1.0, -1.0).normalize()),
            ((0, 0, -1), Vector3::new(-1.0, 1.0, -1.0).normalize()),
        ];

        for (cell_img, bond_dir) in nn_bonds {
            // Projection operator P = r_hat * r_hat^T
            let mut p_proj = [[0.0; 3]; 3];
            let r = [bond_dir.x, bond_dir.y, bond_dir.z];
            for i in 0..3 {
                for j in 0..3 {
                    p_proj[i][j] = r[i] * r[j];
                }
            }

            // Phi = -(alpha * P + beta * (I - P))
            let mut phi = [[0.0; 3]; 3];
            for i in 0..3 {
                for j in 0..3 {
                    let eye = if i == j { 1.0 } else { 0.0 };
                    let val = alpha * p_proj[i][j] + beta * (eye - p_proj[i][j]);
                    phi[i][j] = -val;
                }
            }

            let tensor = Matrix3::new(phi);
            ifc.add_entry_ev_per_angstrom2(0, cell_img, 1, tensor);

            // Symmetrical interaction (1 interacting with 0)
            let inv_cell = (-cell_img.0, -cell_img.1, -cell_img.2);
            ifc.add_entry_ev_per_angstrom2(1, inv_cell, 0, tensor);
        }

        // 2nd nearest neighbors (central force between like atoms on the same sublattice)
        let alpha2 = 0.65;
        let second_nn = [
            (1, 0, 0),
            (-1, 0, 0),
            (0, 1, 0),
            (0, -1, 0),
            (0, 0, 1),
            (0, 0, -1),
        ];

        for cell_img in second_nn {
            let tensor2 = Matrix3::IDENTITY.mul_scalar(-alpha2);
            ifc.add_entry_ev_per_angstrom2(0, cell_img, 0, tensor2);
            ifc.add_entry_ev_per_angstrom2(1, cell_img, 1, tensor2);
        }

        // Enforce acoustic sum rule to guarantee acoustic frequencies == 0 at Gamma
        ifc.enforce_acoustic_sum_rule();

        ifc
    }

    /// Standard FCC Brillouin zone high-symmetry path: Gamma -> X -> W -> L -> Gamma -> K
    pub fn high_symmetry_path(
        reciprocal: &ReciprocalLattice,
        points_per_seg: usize,
    ) -> HighSymmetryPath {
        let waypoints = [
            (Vector3::new(0.0, 0.0, 0.0), "Γ"),
            (Vector3::new(0.5, 0.0, 0.5), "X"),
            (Vector3::new(0.5, 0.25, 0.75), "W"),
            (Vector3::new(0.5, 0.5, 0.5), "L"),
            (Vector3::new(0.0, 0.0, 0.0), "Γ"),
            (Vector3::new(0.375, 0.375, 0.75), "K"),
        ];
        HighSymmetryPath::from_waypoints(reciprocal, &waypoints, points_per_seg)
    }
}
