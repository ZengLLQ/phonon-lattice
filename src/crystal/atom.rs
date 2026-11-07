use crate::math::{Matrix3, Vector3};
use serde::{Deserialize, Serialize};

pub const AMU_TO_KG: f64 = 1.66053906660e-27;

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Atom {
    pub species: String,
    pub mass_amu: f64,
    pub position_fractional: Vector3,
    pub born_charge: Option<Matrix3>,
}

impl Atom {
    pub fn new(species: impl Into<String>, mass_amu: f64, position_fractional: Vector3) -> Self {
        Self {
            species: species.into(),
            mass_amu,
            position_fractional,
            born_charge: None,
        }
    }

    pub fn with_born_charge(mut self, z_star: Matrix3) -> Self {
        self.born_charge = Some(z_star);
        self
    }

    pub fn mass_kg(&self) -> f64 {
        self.mass_amu * AMU_TO_KG
    }
}
