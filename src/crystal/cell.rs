use super::atom::Atom;
use crate::math::{Matrix3, Vector3};
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct UnitCell {
    /// Lattice vectors [a1, a2, a3] in Angstroms (10^-10 m)
    pub a1: Vector3,
    pub a2: Vector3,
    pub a3: Vector3,
    /// Atoms in the primitive cell
    pub atoms: Vec<Atom>,
    /// High-frequency dielectric tensor epsilon_infinity (dimensionless)
    pub dielectric_tensor: Matrix3,
}

impl UnitCell {
    pub fn new(a1: Vector3, a2: Vector3, a3: Vector3, atoms: Vec<Atom>) -> Self {
        Self {
            a1,
            a2,
            a3,
            atoms,
            dielectric_tensor: Matrix3::IDENTITY,
        }
    }

    pub fn with_dielectric_tensor(mut self, eps: Matrix3) -> Self {
        self.dielectric_tensor = eps;
        self
    }

    /// Cell volume in Angstroms^3
    pub fn volume(&self) -> f64 {
        self.a1.dot(&self.a2.cross(&self.a3)).abs()
    }

    /// Cell volume in m^3
    pub fn volume_m3(&self) -> f64 {
        self.volume() * 1e-30
    }

    /// Matrix where columns are lattice vectors A = [a1, a2, a3]
    pub fn lattice_matrix(&self) -> Matrix3 {
        Matrix3::from_cols(self.a1, self.a2, self.a3)
    }

    /// Converts fractional coordinates in [0, 1)^3 to Cartesian coordinates in Angstroms
    pub fn fractional_to_cartesian(&self, frac: Vector3) -> Vector3 {
        self.a1 * frac.x + self.a2 * frac.y + self.a3 * frac.z
    }

    /// Converts Cartesian coordinates in Angstroms to fractional coordinates
    pub fn cartesian_to_fractional(&self, cart: Vector3) -> Option<Vector3> {
        let mat = self.lattice_matrix();
        let inv = mat.inverse()?;
        Some(inv.mul_vec(cart))
    }

    /// Minimum periodic image vector (in fractional coordinates) between frac1 and frac2
    pub fn minimum_image_fractional(frac1: Vector3, frac2: Vector3) -> Vector3 {
        let mut diff = frac1 - frac2;
        diff.x -= diff.x.round();
        diff.y -= diff.y.round();
        diff.z -= diff.z.round();
        diff
    }

    /// Total number of atoms in the unit cell
    pub fn num_atoms(&self) -> usize {
        self.atoms.len()
    }

    /// Total mass of the unit cell in kg
    pub fn total_mass_kg(&self) -> f64 {
        self.atoms.iter().map(|a| a.mass_kg()).sum()
    }

    /// Mass density in kg/m^3
    pub fn mass_density(&self) -> f64 {
        self.total_mass_kg() / self.volume_m3()
    }
}
