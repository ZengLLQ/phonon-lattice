use super::vector3::Vector3;
use serde::{Deserialize, Serialize};
use std::ops::{Add, Mul, Sub};

#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct Matrix3 {
    pub data: [[f64; 3]; 3],
}

impl Matrix3 {
    pub const ZERO: Self = Self {
        data: [[0.0; 3]; 3],
    };

    pub const IDENTITY: Self = Self {
        data: [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
    };

    pub fn new(data: [[f64; 3]; 3]) -> Self {
        Self { data }
    }

    pub fn from_cols(c0: Vector3, c1: Vector3, c2: Vector3) -> Self {
        Self {
            data: [[c0.x, c1.x, c2.x], [c0.y, c1.y, c2.y], [c0.z, c1.z, c2.z]],
        }
    }

    pub fn from_rows(r0: Vector3, r1: Vector3, r2: Vector3) -> Self {
        Self {
            data: [[r0.x, r0.y, r0.z], [r1.x, r1.y, r1.z], [r2.x, r2.y, r2.z]],
        }
    }

    pub fn row(&self, r: usize) -> Vector3 {
        Vector3::new(self.data[r][0], self.data[r][1], self.data[r][2])
    }

    pub fn col(&self, c: usize) -> Vector3 {
        Vector3::new(self.data[0][c], self.data[1][c], self.data[2][c])
    }

    pub fn transpose(&self) -> Self {
        let mut t = [[0.0; 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                t[i][j] = self.data[j][i];
            }
        }
        Self { data: t }
    }

    pub fn determinant(&self) -> f64 {
        let d = &self.data;
        d[0][0] * (d[1][1] * d[2][2] - d[1][2] * d[2][1])
            - d[0][1] * (d[1][0] * d[2][2] - d[1][2] * d[2][0])
            + d[0][2] * (d[1][0] * d[2][1] - d[1][1] * d[2][0])
    }

    pub fn inverse(&self) -> Option<Self> {
        let det = self.determinant();
        if det.abs() < 1e-15 {
            return None;
        }
        let inv_det = 1.0 / det;
        let d = &self.data;

        let inv = [
            [
                (d[1][1] * d[2][2] - d[1][2] * d[2][1]) * inv_det,
                (d[0][2] * d[2][1] - d[0][1] * d[2][2]) * inv_det,
                (d[0][1] * d[1][2] - d[0][2] * d[1][1]) * inv_det,
            ],
            [
                (d[1][2] * d[2][0] - d[1][0] * d[2][2]) * inv_det,
                (d[0][0] * d[2][2] - d[0][2] * d[2][0]) * inv_det,
                (d[0][2] * d[1][0] - d[0][0] * d[1][2]) * inv_det,
            ],
            [
                (d[1][0] * d[2][1] - d[1][1] * d[2][0]) * inv_det,
                (d[0][1] * d[2][0] - d[0][0] * d[2][1]) * inv_det,
                (d[0][0] * d[1][1] - d[0][1] * d[1][0]) * inv_det,
            ],
        ];

        Some(Self { data: inv })
    }

    pub fn mul_vec(&self, v: Vector3) -> Vector3 {
        Vector3::new(
            self.data[0][0] * v.x + self.data[0][1] * v.y + self.data[0][2] * v.z,
            self.data[1][0] * v.x + self.data[1][1] * v.y + self.data[1][2] * v.z,
            self.data[2][0] * v.x + self.data[2][1] * v.y + self.data[2][2] * v.z,
        )
    }

    pub fn mul_scalar(&self, s: f64) -> Self {
        let mut res = [[0.0; 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                res[i][j] = self.data[i][j] * s;
            }
        }
        Self { data: res }
    }

    pub fn trace(&self) -> f64 {
        self.data[0][0] + self.data[1][1] + self.data[2][2]
    }
}

impl Add for Matrix3 {
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        let mut res = [[0.0; 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                res[i][j] = self.data[i][j] + rhs.data[i][j];
            }
        }
        Self { data: res }
    }
}

impl Sub for Matrix3 {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self {
        let mut res = [[0.0; 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                res[i][j] = self.data[i][j] - rhs.data[i][j];
            }
        }
        Self { data: res }
    }
}

impl Mul for Matrix3 {
    type Output = Self;
    fn mul(self, rhs: Self) -> Self {
        let mut res = [[0.0; 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                let mut sum = 0.0;
                for k in 0..3 {
                    sum += self.data[i][k] * rhs.data[k][j];
                }
                res[i][j] = sum;
            }
        }
        Self { data: res }
    }
}
