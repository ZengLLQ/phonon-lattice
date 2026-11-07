use super::dynamical_matrix::DynamicalMatrixBuilder;
use crate::crystal::HighSymmetryPath;
use crate::math::{eigensystem_hermitian, Vector3};
use serde::{Deserialize, Serialize};
use std::f64::consts::PI;

pub const RAD_S_TO_THZ: f64 = 1.0 / (2.0 * PI * 1e12);
pub const THZ_TO_CM_MINUS_1: f64 = 33.35641;
pub const THZ_TO_MEV: f64 = 4.135667696;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DispersionPoint {
    pub fractional_q: Vector3,
    pub cartesian_q: Vector3,
    pub cumulative_distance: f64,
    /// Phonon frequencies for each branch s in THz
    pub frequencies_thz: Vec<f64>,
    /// Phonon frequencies in rad/s
    pub frequencies_rad_s: Vec<f64>,
    /// Group velocities in m/s for each branch
    pub group_velocities_m_s: Vec<Vector3>,
    /// True if an imaginary / dynamically soft mode was detected
    pub has_imaginary_modes: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DispersionBand {
    pub branch_index: usize,
    pub label: String,
    pub frequencies_thz: Vec<f64>,
    pub max_frequency_thz: f64,
    pub min_frequency_thz: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PhononDispersion {
    pub num_branches: usize,
    pub path_points: Vec<DispersionPoint>,
    pub special_point_indices: Vec<(usize, String)>,
}

impl PhononDispersion {
    pub fn calculate(builder: &DynamicalMatrixBuilder, path: &HighSymmetryPath) -> Self {
        let n_atoms = builder.cell.num_atoms();
        let num_branches = 3 * n_atoms;

        let mut path_points = Vec::with_capacity(path.points.len());

        for pt in &path.points {
            let d_mat = builder.build_at(pt.cartesian);
            let eigen = eigensystem_hermitian(&d_mat, 1e-13, 200);

            let mut freqs_rad_s = Vec::with_capacity(num_branches);
            let mut freqs_thz = Vec::with_capacity(num_branches);
            let mut has_soft = false;

            for &val in &eigen.eigenvalues {
                if val >= 0.0 {
                    let omega = val.sqrt();
                    freqs_rad_s.push(omega);
                    freqs_thz.push(omega * RAD_S_TO_THZ);
                } else {
                    // Imaginary soft mode (dynamically unstable)
                    let omega = (-val).sqrt();
                    freqs_rad_s.push(-omega);
                    freqs_thz.push(-omega * RAD_S_TO_THZ);
                    has_soft = true;
                }
            }

            let group_vels =
                builder.compute_group_velocities(pt.cartesian, &freqs_rad_s, &eigen.eigenvectors);

            path_points.push(DispersionPoint {
                fractional_q: pt.fractional,
                cartesian_q: pt.cartesian,
                cumulative_distance: pt.cumulative_distance,
                frequencies_thz: freqs_thz,
                frequencies_rad_s: freqs_rad_s,
                group_velocities_m_s: group_vels,
                has_imaginary_modes: has_soft,
            });
        }

        Self {
            num_branches,
            path_points,
            special_point_indices: path.special_point_indices.clone(),
        }
    }

    /// Maximum phonon frequency across all branches in THz
    pub fn max_frequency_thz(&self) -> f64 {
        self.path_points
            .iter()
            .flat_map(|p| p.frequencies_thz.iter())
            .cloned()
            .fold(0.0, f64::max)
    }

    /// Extract a single branch across the dispersion path
    pub fn extract_band(&self, branch_idx: usize) -> DispersionBand {
        let freqs: Vec<f64> = self
            .path_points
            .iter()
            .map(|p| p.frequencies_thz[branch_idx])
            .collect();
        let max_f = freqs.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let min_f = freqs.iter().cloned().fold(f64::INFINITY, f64::min);
        let label = if branch_idx < 3 {
            format!("Acoustic #{}", branch_idx + 1)
        } else {
            format!("Optical #{}", branch_idx - 2)
        };
        DispersionBand {
            branch_index: branch_idx,
            label,
            frequencies_thz: freqs,
            max_frequency_thz: max_f,
            min_frequency_thz: min_f,
        }
    }

    /// Renders an ASCII text-based dispersion plot
    pub fn render_ascii_plot(&self, height: usize, width: usize) -> String {
        let max_f = self.max_frequency_thz().max(1.0);
        let min_f = 0.0;
        let n_pts = self.path_points.len();
        if n_pts == 0 {
            return String::new();
        }

        let mut canvas = vec![vec![' '; width]; height];

        for col in 0..width {
            let pt_idx = (col * (n_pts - 1)) / (width - 1);
            let pt = &self.path_points[pt_idx];

            for &f in &pt.frequencies_thz {
                if f >= min_f && f <= max_f {
                    let row = ((1.0 - (f - min_f) / (max_f - min_f)) * (height as f64 - 1.0))
                        .round() as usize;
                    let r = row.min(height - 1);
                    canvas[r][col] = '*';
                }
            }
        }

        // Draw vertical guide lines for special points
        for &(idx, _) in &self.special_point_indices {
            let col = (idx * (width - 1)) / (n_pts - 1);
            for r in 0..height {
                if canvas[r][col] == ' ' {
                    canvas[r][col] = '|';
                }
            }
        }

        let mut out = String::new();
        for r in 0..height {
            let freq = max_f * (1.0 - r as f64 / (height as f64 - 1.0));
            out.push_str(&format!("{:5.1} THz |", freq));
            for c in 0..width {
                out.push(canvas[r][c]);
            }
            out.push('\n');
        }

        out.push_str("          +");
        for _ in 0..width {
            out.push('-');
        }
        out.push('\n');

        // Bottom label line
        out.push_str("           ");
        let mut label_line = vec![' '; width];
        for &(idx, ref name) in &self.special_point_indices {
            let col = (idx * (width - 1)) / (n_pts - 1);
            let chars: Vec<char> = name.chars().collect();
            for (offset, &ch) in chars.iter().enumerate() {
                if col + offset < width {
                    label_line[col + offset] = ch;
                }
            }
        }
        for ch in label_line {
            out.push(ch);
        }
        out.push('\n');

        out
    }
}
