use num_complex::Complex64;

#[derive(Debug, Clone)]
pub struct EigenResult {
    /// Eigenvalues sorted in ascending order
    pub eigenvalues: Vec<f64>,
    /// Column eigenvectors corresponding to each eigenvalue.
    /// eigenvectors[j] is the j-th eigenvector (length N).
    pub eigenvectors: Vec<Vec<Complex64>>,
}

/// Solves the eigenvalue problem for an N x N complex Hermitian matrix A:
/// A * v_j = lambda_j * v_j
/// where A[i][j] == conj(A[j][i]).
///
/// Uses the complex Jacobi cyclic rotation algorithm with quadratic convergence.
pub fn eigensystem_hermitian(matrix: &[Vec<Complex64>], tol: f64, max_iter: usize) -> EigenResult {
    let n = matrix.len();
    assert!(n > 0, "Matrix dimension must be positive");
    for row in matrix {
        assert_eq!(row.len(), n, "Matrix must be square");
    }

    // Work matrix initialized to A
    let mut a = matrix.to_vec();

    // Enforce exact Hermiticity on input
    for i in 0..n {
        a[i][i] = Complex64::new(a[i][i].re, 0.0);
        for j in (i + 1)..n {
            let avg = (a[i][j] + a[j][i].conj()) * 0.5;
            a[i][j] = avg;
            a[j][i] = avg.conj();
        }
    }

    // Accumulator for unitary eigenvectors, initialized to Identity matrix
    let mut v = vec![vec![Complex64::new(0.0, 0.0); n]; n];
    for i in 0..n {
        v[i][i] = Complex64::new(1.0, 0.0);
    }

    for _iter in 0..max_iter {
        // Calculate sum of squares of off-diagonal elements
        let mut off_diag_sum = 0.0;
        for i in 0..n {
            for j in (i + 1)..n {
                off_diag_sum += a[i][j].norm_sqr();
            }
        }

        if off_diag_sum.sqrt() < tol {
            break;
        }

        // Cyclic sweep over all off-diagonal pairs (p, q)
        for p in 0..(n - 1) {
            for q in (p + 1)..n {
                let apq = a[p][q];
                let mag = apq.norm();
                if mag < 1e-18 {
                    continue;
                }

                // Phase angle phi of A_pq = |A_pq| * exp(i * phi)
                let phase = apq / mag;

                // 2x2 submatrix:
                // [ a_pp,   |a_pq| ]
                // [ |a_pq|, a_qq   ]
                let app = a[p][p].re;
                let aqq = a[q][q].re;
                let tau = (aqq - app) / (2.0 * mag);

                let t = if tau >= 0.0 {
                    1.0 / (tau + (1.0 + tau * tau).sqrt())
                } else {
                    -1.0 / (-tau + (1.0 + tau * tau).sqrt())
                };

                let c = 1.0 / (1.0 + t * t).sqrt();
                let s = t * c;

                // Jacobi rotation parameters:
                // U_pp = c, U_pq = s * phase, U_qp = -s * conj(phase), U_qq = c
                let s_phase = phase * s;
                let s_phase_conj = phase.conj() * s;

                // Update diagonal elements
                a[p][p] = Complex64::new(app - t * mag, 0.0);
                a[q][q] = Complex64::new(aqq + t * mag, 0.0);
                a[p][q] = Complex64::new(0.0, 0.0);
                a[q][p] = Complex64::new(0.0, 0.0);

                // Update elements in rows/columns other than p and q
                for r in 0..n {
                    if r != p && r != q {
                        let arp = a[r][p];
                        let arq = a[r][q];

                        // A' = U^\dagger * A * U
                        let new_arp = arp * c - arq * s_phase_conj;
                        let new_arq = arq * c + arp * s_phase;

                        a[r][p] = new_arp;
                        a[p][r] = new_arp.conj();
                        a[r][q] = new_arq;
                        a[q][r] = new_arq.conj();
                    }
                }

                // Update accumulated eigenvector matrix V: V' = V * U
                for r in 0..n {
                    let vrp = v[r][p];
                    let vrq = v[r][q];
                    v[r][p] = vrp * c - vrq * s_phase_conj;
                    v[r][q] = vrq * c + vrp * s_phase;
                }
            }
        }
    }

    // Extract eigenvalues and pack eigenvectors
    let mut pairs: Vec<(f64, Vec<Complex64>)> = (0..n)
        .map(|j| {
            let lambda = a[j][j].re;
            let vec_j = (0..n).map(|i| v[i][j]).collect();
            (lambda, vec_j)
        })
        .collect();

    // Sort by eigenvalue ascending
    pairs.sort_by(|(l1, _), (l2, _)| l1.partial_cmp(l2).unwrap_or(std::cmp::Ordering::Equal));

    let eigenvalues = pairs.iter().map(|(l, _)| *l).collect();
    let eigenvectors = pairs.into_iter().map(|(_, vec)| vec).collect();

    EigenResult {
        eigenvalues,
        eigenvectors,
    }
}
