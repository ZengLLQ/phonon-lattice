use crate::crystal::{Atom, HighSymmetryPath, ReciprocalLattice, UnitCell};
use crate::dynamics::force_constants::InteratomicForceConstants;
use crate::math::{Matrix3, Vector3};

pub struct GrapheneMaterial;

impl GrapheneMaterial {
    pub const LATTICE_CONSTANT_A: f64 = 2.46; // Angstroms
    pub const VACUUM_C: f64 = 20.0; // Angstroms
    pub const MASS_CARBON_AMU: f64 = 12.011;

    pub fn unit_cell() -> UnitCell {
        let a = Self::LATTICE_CONSTANT_A;
        let a1 = Vector3::new(a * 3.0f64.sqrt() * 0.5, -a * 0.5, 0.0);
        let a2 = Vector3::new(a * 3.0f64.sqrt() * 0.5, a * 0.5, 0.0);
        let a3 = Vector3::new(0.0, 0.0, Self::VACUUM_C);

        let atoms = vec![
            Atom::new(
                "C",
                Self::MASS_CARBON_AMU,
                Vector3::new(1.0 / 3.0, 2.0 / 3.0, 0.0),
            ),
            Atom::new(
                "C",
                Self::MASS_CARBON_AMU,
                Vector3::new(2.0 / 3.0, 1.0 / 3.0, 0.0),
            ),
        ];

        UnitCell::new(a1, a2, a3, atoms)
    }

    pub fn force_constants() -> InteratomicForceConstants {
        let mut ifc = InteratomicForceConstants::new(2);

        // Carbon-Carbon bond parameters (eV/Angstrom^2)
        // alpha_r: radial bond-stretching (~25.8 eV/A^2)
        // beta_ti: in-plane transverse bond-bending (~9.2 eV/A^2)
        // gamma_to: out-of-plane flexural bond-bending (~5.8 eV/A^2)
        let alpha_r = 25.8;
        let beta_ti = 9.2;
        let gamma_to = 5.8;

        // 3 nearest neighbors in honeycomb lattice for atom 0
        let nn_bonds = [
            ((0, 0, 0), Vector3::new(1.0, 0.0, 0.0)),
            ((-1, 0, 0), Vector3::new(-0.5, 3.0f64.sqrt() * 0.5, 0.0)),
            ((0, -1, 0), Vector3::new(-0.5, -3.0f64.sqrt() * 0.5, 0.0)),
        ];

        for (cell_img, bond_dir) in nn_bonds {
            let r_in_plane = bond_dir.normalize();
            let t_in_plane = Vector3::new(-r_in_plane.y, r_in_plane.x, 0.0);
            let z_out = Vector3::new(0.0, 0.0, 1.0);

            let mut phi = [[0.0; 3]; 3];
            let r = [r_in_plane.x, r_in_plane.y, r_in_plane.z];
            let t = [t_in_plane.x, t_in_plane.y, t_in_plane.z];
            let z = [z_out.x, z_out.y, z_out.z];

            for i in 0..3 {
                for j in 0..3 {
                    phi[i][j] =
                        -(alpha_r * r[i] * r[j] + beta_ti * t[i] * t[j] + gamma_to * z[i] * z[j]);
                }
            }

            let tensor = Matrix3::new(phi);
            ifc.add_entry_ev_per_angstrom2(0, cell_img, 1, tensor);

            let inv_cell = (-cell_img.0, -cell_img.1, -cell_img.2);
            ifc.add_entry_ev_per_angstrom2(1, inv_cell, 0, tensor);
        }

        // 2nd nearest neighbor in-plane interaction
        let alpha2 = 1.4;
        let second_nn = [
            (1, 0, 0),
            (-1, 0, 0),
            (0, 1, 0),
            (0, -1, 0),
            (1, -1, 0),
            (-1, 1, 0),
        ];

        for cell_img in second_nn {
            let t2 = Matrix3::new([
                [-alpha2, 0.0, 0.0],
                [0.0, -alpha2, 0.0],
                [0.0, 0.0, -alpha2 * 0.3],
            ]);
            ifc.add_entry_ev_per_angstrom2(0, cell_img, 0, t2);
            ifc.add_entry_ev_per_angstrom2(1, cell_img, 1, t2);
        }

        ifc.enforce_acoustic_sum_rule();
        ifc
    }

    /// Hexagonal Brillouin zone path: Gamma -> M -> K -> Gamma
    pub fn high_symmetry_path(
        reciprocal: &ReciprocalLattice,
        points_per_seg: usize,
    ) -> HighSymmetryPath {
        let waypoints = [
            (Vector3::new(0.0, 0.0, 0.0), "Γ"),
            (Vector3::new(0.5, 0.0, 0.0), "M"),
            (Vector3::new(1.0 / 3.0, 1.0 / 3.0, 0.0), "K"),
            (Vector3::new(0.0, 0.0, 0.0), "Γ"),
        ];
        HighSymmetryPath::from_waypoints(reciprocal, &waypoints, points_per_seg)
    }
}
