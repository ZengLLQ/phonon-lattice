use super::cell::UnitCell;
use crate::math::{Matrix3, Vector3};
use serde::{Deserialize, Serialize};
use std::f64::consts::PI;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ReciprocalLattice {
    /// Reciprocal vectors b1, b2, b3 in Angstroms^-1
    pub b1: Vector3,
    pub b2: Vector3,
    pub b3: Vector3,
}

impl ReciprocalLattice {
    pub fn from_unit_cell(cell: &UnitCell) -> Self {
        let vol = cell.volume();
        let factor = 2.0 * PI / vol;
        Self {
            b1: cell.a2.cross(&cell.a3) * factor,
            b2: cell.a3.cross(&cell.a1) * factor,
            b3: cell.a1.cross(&cell.a2) * factor,
        }
    }

    /// Matrix whose columns are reciprocal basis vectors B = [b1, b2, b3]
    pub fn reciprocal_matrix(&self) -> Matrix3 {
        Matrix3::from_cols(self.b1, self.b2, self.b3)
    }

    /// Converts fractional reciprocal coordinates (k1, k2, k3) to Cartesian q in Angstroms^-1
    pub fn fractional_to_cartesian(&self, frac: Vector3) -> Vector3 {
        self.b1 * frac.x + self.b2 * frac.y + self.b3 * frac.z
    }

    /// Converts Cartesian q in Angstroms^-1 to fractional reciprocal coordinates
    pub fn cartesian_to_fractional(&self, cart: Vector3) -> Option<Vector3> {
        let mat = self.reciprocal_matrix();
        let inv = mat.inverse()?;
        Some(inv.mul_vec(cart))
    }

    /// Generates a uniform Monkhorst-Pack grid of fractional q-points centered in [-0.5, 0.5)^3
    pub fn generate_monkhorst_pack_grid(&self, n1: usize, n2: usize, n3: usize) -> Vec<Vector3> {
        let mut q_points = Vec::with_capacity(n1 * n2 * n3);
        for i in 0..n1 {
            let u1 = (2 * i as i32 - n1 as i32 + 1) as f64 / (2.0 * n1 as f64);
            for j in 0..n2 {
                let u2 = (2 * j as i32 - n2 as i32 + 1) as f64 / (2.0 * n2 as f64);
                for k in 0..n3 {
                    let u3 = (2 * k as i32 - n3 as i32 + 1) as f64 / (2.0 * n3 as f64);
                    q_points.push(Vector3::new(u1, u2, u3));
                }
            }
        }
        q_points
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PathPoint {
    pub fractional: Vector3,
    pub cartesian: Vector3,
    pub label: Option<String>,
    pub cumulative_distance: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct HighSymmetryPath {
    pub points: Vec<PathPoint>,
    pub special_point_indices: Vec<(usize, String)>,
}

impl HighSymmetryPath {
    /// Constructs a continuous path across a list of named high-symmetry waypoints.
    /// Each segment between waypoints is sampled with `points_per_segment` equidistant points.
    pub fn from_waypoints(
        reciprocal: &ReciprocalLattice,
        waypoints: &[(Vector3, &str)],
        points_per_segment: usize,
    ) -> Self {
        assert!(
            waypoints.len() >= 2,
            "Path must contain at least 2 waypoints"
        );
        let mut points = Vec::new();
        let mut special_indices = Vec::new();
        let mut cum_dist = 0.0;

        for seg_idx in 0..(waypoints.len() - 1) {
            let (start_frac, start_label) = waypoints[seg_idx];
            let (end_frac, end_label) = waypoints[seg_idx + 1];

            let start_cart = reciprocal.fractional_to_cartesian(start_frac);
            let end_cart = reciprocal.fractional_to_cartesian(end_frac);
            let seg_vec = end_cart - start_cart;
            let seg_len = seg_vec.norm();

            // Record start of segment
            if seg_idx == 0 {
                points.push(PathPoint {
                    fractional: start_frac,
                    cartesian: start_cart,
                    label: Some(start_label.to_string()),
                    cumulative_distance: cum_dist,
                });
                special_indices.push((0, start_label.to_string()));
            }

            let n_steps = points_per_segment.max(2);
            for step in 1..=n_steps {
                let t = step as f64 / n_steps as f64;
                let cur_frac = start_frac * (1.0 - t) + end_frac * t;
                let cur_cart = start_cart * (1.0 - t) + end_cart * t;
                let cur_dist = cum_dist + seg_len * t;

                let is_end = step == n_steps;
                let label = if is_end {
                    Some(end_label.to_string())
                } else {
                    None
                };

                let idx = points.len();
                points.push(PathPoint {
                    fractional: cur_frac,
                    cartesian: cur_cart,
                    label,
                    cumulative_distance: cur_dist,
                });

                if is_end {
                    special_indices.push((idx, end_label.to_string()));
                }
            }

            cum_dist += seg_len;
        }

        Self {
            points,
            special_point_indices: special_indices,
        }
    }
}
