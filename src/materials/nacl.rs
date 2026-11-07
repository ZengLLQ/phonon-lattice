use crate::crystal::{Atom, HighSymmetryPath, ReciprocalLattice, UnitCell};
use crate::dynamics::force_constants::InteratomicForceConstants;
use crate::math::{Matrix3, Vector3};

pub struct NaClMaterial;

impl NaClMaterial {
    pub const LATTICE_CONSTANT: f64 = 5.64; // Angstroms
    pub const MASS_NA_AMU: f64 = 22.990;
    pub const MASS_CL_AMU: f64 = 35.453;
    pub const BORN_CHARGE_Z_STAR: f64 = 1.09;
    pub const DIELECTRIC_EPS_INF: f64 = 2.34;

    pub fn unit_cell() -> UnitCell {
        let a = Self::LATTICE_CONSTANT * 0.5;
        let a1 = Vector3::new(0.0, a, a);
        let a2 = Vector3::new(a, 0.0, a);
        let a3 = Vector3::new(a, a, 0.0);

        let z_na = Matrix3::IDENTITY.mul_scalar(Self::BORN_CHARGE_Z_STAR);
        let z_cl = Matrix3::IDENTITY.mul_scalar(-Self::BORN_CHARGE_Z_STAR);

        let atoms = vec![
            Atom::new("Na", Self::MASS_NA_AMU, Vector3::new(0.0, 0.0, 0.0)).with_born_charge(z_na),
            Atom::new("Cl", Self::MASS_CL_AMU, Vector3::new(0.5, 0.5, 0.5)).with_born_charge(z_cl),
        ];

        let eps = Matrix3::IDENTITY.mul_scalar(Self::DIELECTRIC_EPS_INF);
        UnitCell::new(a1, a2, a3, atoms).with_dielectric_tensor(eps)
    }

    pub fn force_constants() -> InteratomicForceConstants {
        let mut ifc = InteratomicForceConstants::new(2);

        // 6 nearest neighbors in rocksalt octahedron
        let alpha = 2.85; // eV / Angstrom^2
        let beta = 0.55;

        let nn_bonds = [
            ((0, 0, 0), Vector3::new(1.0, 0.0, 0.0)),
            ((-1, 0, 0), Vector3::new(-1.0, 0.0, 0.0)),
            ((0, 0, 0), Vector3::new(0.0, 1.0, 0.0)),
            ((0, -1, 0), Vector3::new(0.0, -1.0, 0.0)),
            ((0, 0, 0), Vector3::new(0.0, 0.0, 1.0)),
            ((0, 0, -1), Vector3::new(0.0, 0.0, -1.0)),
        ];

        for (cell_img, bond_dir) in nn_bonds {
            let mut p_proj = [[0.0; 3]; 3];
            let r = [bond_dir.x, bond_dir.y, bond_dir.z];
            for i in 0..3 {
                for j in 0..3 {
                    p_proj[i][j] = r[i] * r[j];
                }
            }

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

            let inv_cell = (-cell_img.0, -cell_img.1, -cell_img.2);
            ifc.add_entry_ev_per_angstrom2(1, inv_cell, 0, tensor);
        }

        ifc.enforce_acoustic_sum_rule();
        ifc
    }

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
        ];
        HighSymmetryPath::from_waypoints(reciprocal, &waypoints, points_per_seg)
    }
}
