use crate::crystal::{Atom, HighSymmetryPath, ReciprocalLattice, UnitCell};
use crate::dynamics::force_constants::InteratomicForceConstants;
use crate::math::{Matrix3, Vector3};

pub struct GaAsMaterial;

impl GaAsMaterial {
    pub const LATTICE_CONSTANT: f64 = 5.653; // Angstroms
    pub const MASS_GA_AMU: f64 = 69.723;
    pub const MASS_AS_AMU: f64 = 74.922;
    pub const BORN_CHARGE_Z_STAR: f64 = 2.18; // elementary charges
    pub const DIELECTRIC_EPS_INF: f64 = 10.9;

    pub fn unit_cell() -> UnitCell {
        let a = Self::LATTICE_CONSTANT * 0.5;
        let a1 = Vector3::new(0.0, a, a);
        let a2 = Vector3::new(a, 0.0, a);
        let a3 = Vector3::new(a, a, 0.0);

        let z_ga = Matrix3::IDENTITY.mul_scalar(Self::BORN_CHARGE_Z_STAR);
        let z_as = Matrix3::IDENTITY.mul_scalar(-Self::BORN_CHARGE_Z_STAR);

        let atoms = vec![
            Atom::new("Ga", Self::MASS_GA_AMU, Vector3::new(0.0, 0.0, 0.0)).with_born_charge(z_ga),
            Atom::new("As", Self::MASS_AS_AMU, Vector3::new(0.25, 0.25, 0.25))
                .with_born_charge(z_as),
        ];

        let eps = Matrix3::IDENTITY.mul_scalar(Self::DIELECTRIC_EPS_INF);
        UnitCell::new(a1, a2, a3, atoms).with_dielectric_tensor(eps)
    }

    pub fn force_constants() -> InteratomicForceConstants {
        let mut ifc = InteratomicForceConstants::new(2);

        let alpha = 7.15;
        let beta = 2.45;

        let nn_bonds = [
            ((0, 0, 0), Vector3::new(1.0, 1.0, 1.0).normalize()),
            ((-1, 0, 0), Vector3::new(-1.0, -1.0, 1.0).normalize()),
            ((0, -1, 0), Vector3::new(1.0, -1.0, -1.0).normalize()),
            ((0, 0, -1), Vector3::new(-1.0, 1.0, -1.0).normalize()),
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

        let alpha2 = 0.42;
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
