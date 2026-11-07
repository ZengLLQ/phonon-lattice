use num_complex::Complex64;
use phonon_lattice::math::eigensystem_hermitian;

#[test]
fn test_hermitian_pauli_y_matrix() {
    // Pauli sigma_y = [[0, -i], [i, 0]]
    // Eigenvalues should be -1.0 and +1.0
    let sigma_y = vec![
        vec![Complex64::new(0.0, 0.0), Complex64::new(0.0, -1.0)],
        vec![Complex64::new(0.0, 1.0), Complex64::new(0.0, 0.0)],
    ];

    let result = eigensystem_hermitian(&sigma_y, 1e-14, 50);

    assert_eq!(result.eigenvalues.len(), 2);
    assert!((result.eigenvalues[0] - (-1.0)).abs() < 1e-12);
    assert!((result.eigenvalues[1] - 1.0).abs() < 1e-12);

    // Verify eigenvector orthogonality: <v_0 | v_1> == 0
    let dot = result.eigenvectors[0][0].conj() * result.eigenvectors[1][0]
        + result.eigenvectors[0][1].conj() * result.eigenvectors[1][1];
    assert!(dot.norm() < 1e-12);
}

#[test]
fn test_hermitian_3x3_arbitrary_matrix() {
    let mat = vec![
        vec![
            Complex64::new(5.0, 0.0),
            Complex64::new(2.0, 1.0),
            Complex64::new(0.0, -3.0),
        ],
        vec![
            Complex64::new(2.0, -1.0),
            Complex64::new(4.0, 0.0),
            Complex64::new(1.0, 2.0),
        ],
        vec![
            Complex64::new(0.0, 3.0),
            Complex64::new(1.0, -2.0),
            Complex64::new(7.0, 0.0),
        ],
    ];

    let result = eigensystem_hermitian(&mat, 1e-14, 100);

    // Check trace conservation: sum(lambda) == tr(A) = 5 + 4 + 7 = 16.0
    let trace_sum: f64 = result.eigenvalues.iter().sum();
    assert!((trace_sum - 16.0).abs() < 1e-10);

    // Check eigenvalue equation: A * v_j = lambda_j * v_j
    for j in 0..3 {
        let lambda = result.eigenvalues[j];
        let v = &result.eigenvectors[j];

        for row in 0..3 {
            let mut lhs = Complex64::new(0.0, 0.0);
            for col in 0..3 {
                lhs += mat[row][col] * v[col];
            }
            let rhs = v[row] * lambda;
            assert!(
                (lhs - rhs).norm() < 1e-10,
                "Eigenvalue equation failed for mode {}, row {}: lhs={:?}, rhs={:?}",
                j,
                row,
                lhs,
                rhs
            );
        }
    }
}
