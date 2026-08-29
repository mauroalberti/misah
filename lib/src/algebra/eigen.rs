//! Eigen-decomposition of a symmetric 3x3 matrix, by cyclic Jacobi rotations.
//!
//! Exists because the best-fit plane needs the direction of least scatter in a
//! cloud of points, and that is the eigenvector of the smallest eigenvalue of
//! the cloud's scatter matrix. The reference implementation this is ported
//! alongside reached the same answer through LAPACK's `dgesvd` on the m-by-3
//! matrix of centred coordinates; the right singular vectors of that matrix
//! are the eigenvectors of its scatter matrix, and the singular values the
//! square roots of the eigenvalues, so the two routes agree by construction.
//!
//! Taking the second route costs one thing and saves another. It squares the
//! condition number: forming `X^T X` before decomposing it loses about half
//! the available significant digits, where an SVD of `X` would not. What it
//! saves is LAPACK and BLAS, which is the whole numerical-library dependency
//! this crate does not otherwise have, for a matrix small enough to decompose
//! exactly in a page of code. For plane fitting the trade is not close: the
//! coordinates are centred first, so the entries are metres of local relief
//! rather than full projected coordinates, and an attitude is wanted to a
//! hundredth of a degree, not to the fifteenth digit.
//!
//! Jacobi rather than a closed form. A symmetric 3x3 has an analytic
//! eigen-solution through the characteristic cubic, and it is faster; it is
//! also the one that loses accuracy exactly where this is used most, on nearly
//! degenerate matrices, where two eigenvalues are close and the cubic's roots
//! are ill-conditioned even though the eigenvectors are not.

/// Eigenvalues in descending order, and their eigenvectors as the columns of
/// the returned matrix: `eigenvectors[row][k]` is component `row` of the
/// eigenvector for `eigenvalues[k]`.
///
/// The matrix is assumed symmetric and only its upper triangle is read.
/// Eigenvectors are orthonormal, and each is determined only up to sign --
/// which is not a shortcoming here, since the direction of a plane's normal
/// carries no up or down either.
pub fn symmetric_eigen(matrix: [[f64; 3]; 3]) -> ([f64; 3], [[f64; 3]; 3]) {

    let mut a = matrix;
    let mut v = [
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
    ];

    // Symmetry is assumed, not checked: mirror the upper triangle down so that
    // a caller who filled in only half is decomposed rather than silently
    // given the answer to a different matrix.
    a[1][0] = a[0][1];
    a[2][0] = a[0][2];
    a[2][1] = a[1][2];

    // Each sweep annihilates the largest off-diagonal entry. Convergence is
    // quadratic once the matrix is nearly diagonal, so a 3x3 is done in well
    // under a dozen; the cap is there to bound the loop, not because it is
    // expected to be reached.
    for _ in 0..MAX_ROTATIONS {

        let (p, q, largest) = largest_off_diagonal(&a);

        // Scaled by the diagonal, so the test means "small compared with this
        // matrix" rather than "small in whatever units the caller used".
        let scale = a[0][0].abs() + a[1][1].abs() + a[2][2].abs();
        if largest <= f64::EPSILON * scale || largest == 0.0 {
            break;
        }

        // The rotation angle that zeroes a[p][q], in the stable form: taking
        // the root with the smaller magnitude keeps the rotation the lesser of
        // the two that would do it, which is what makes repeated sweeps
        // converge rather than wander.
        let theta = (a[q][q] - a[p][p]) / (2.0 * a[p][q]);
        let t = theta.signum() / (theta.abs() + (theta * theta + 1.0).sqrt());
        let c = 1.0 / (t * t + 1.0).sqrt();
        let s = t * c;

        rotate(&mut a, &mut v, p, q, c, s);
    }

    let mut eigenvalues = [a[0][0], a[1][1], a[2][2]];

    // Descending, carrying each eigenvector along with its own eigenvalue.
    // Three elements, so the sort is written out: two passes of the largest
    // remaining to the front.
    for i in 0..2 {
        let mut largest = i;
        for j in (i + 1)..3 {
            if eigenvalues[j] > eigenvalues[largest] {
                largest = j;
            }
        }
        if largest != i {
            eigenvalues.swap(i, largest);
            for row in v.iter_mut() {
                row.swap(i, largest);
            }
        }
    }

    (eigenvalues, v)
}

/// Enough rotations for a 3x3 to converge many times over, and a bound on the
/// loop if a matrix carrying NaN makes the tolerance test unsatisfiable.
const MAX_ROTATIONS: usize = 50;

/// The off-diagonal entry of largest magnitude, as `(row, column, magnitude)`.
fn largest_off_diagonal(a: &[[f64; 3]; 3]) -> (usize, usize, f64) {

    let mut best = (0usize, 1usize, a[0][1].abs());

    for &(p, q) in &[(0usize, 2usize), (1, 2)] {
        if a[p][q].abs() > best.2 {
            best = (p, q, a[p][q].abs());
        }
    }

    best
}

/// Apply one Jacobi rotation in the (p, q) plane to `a`, accumulating it in `v`.
fn rotate(a: &mut [[f64; 3]; 3], v: &mut [[f64; 3]; 3], p: usize, q: usize, c: f64, s: f64) {

    let (app, aqq, apq) = (a[p][p], a[q][q], a[p][q]);

    // The two diagonal entries take the whole of the annihilated off-diagonal
    // between them, which is why the trace is preserved exactly.
    a[p][p] = c * c * app - 2.0 * s * c * apq + s * s * aqq;
    a[q][q] = s * s * app + 2.0 * s * c * apq + c * c * aqq;
    a[p][q] = 0.0;
    a[q][p] = 0.0;

    // The remaining row and column: the one index that is neither p nor q.
    let r = 3 - p - q;
    let (arp, arq) = (a[r][p], a[r][q]);
    a[r][p] = c * arp - s * arq;
    a[p][r] = a[r][p];
    a[r][q] = s * arp + c * arq;
    a[q][r] = a[r][q];

    for row in v.iter_mut() {
        let (rp, rq) = (row[p], row[q]);
        row[p] = c * rp - s * rq;
        row[q] = s * rp + c * rq;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// `V . diag(values) . V^T`, to put back together what was taken apart.
    fn recompose(values: [f64; 3], vectors: [[f64; 3]; 3]) -> [[f64; 3]; 3] {

        let mut out = [[0.0; 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                out[i][j] = (0..3).map(|k| vectors[i][k] * values[k] * vectors[j][k]).sum();
            }
        }
        out
    }

    /// Inclusive, so that a tolerance of zero asserts exact equality rather
    /// than asserting nothing can pass.
    fn assert_close(a: [[f64; 3]; 3], b: [[f64; 3]; 3], tolerance: f64) {
        for i in 0..3 {
            for j in 0..3 {
                assert!(
                    (a[i][j] - b[i][j]).abs() <= tolerance,
                    "[{i}][{j}]: {} against {}",
                    a[i][j],
                    b[i][j]
                );
            }
        }
    }

    #[test]
    fn a_diagonal_matrix_is_already_decomposed() {
        let (values, vectors) = symmetric_eigen([[3.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 2.0]]);

        assert_eq!(values, [3.0, 2.0, 1.0]);
        // Sorted descending, so the eigenvector for 2 is the third axis.
        assert!((vectors[2][1].abs() - 1.0).abs() < 1e-12);
    }

    #[test]
    fn the_identity_decomposes_to_itself() {
        let (values, vectors) = symmetric_eigen([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]);

        assert_eq!(values, [1.0, 1.0, 1.0]);
        // Every direction is an eigenvector, so only orthonormality is
        // meaningful; the vectors themselves are whatever the sweep left.
        assert_close(recompose(values, vectors), [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], 1e-12);
    }

    #[test]
    fn the_zero_matrix_does_not_spin() {
        // Every off-diagonal is already zero, so the loop must exit on its
        // first test rather than on the rotation cap.
        let (values, vectors) = symmetric_eigen([[0.0; 3]; 3]);

        assert_eq!(values, [0.0, 0.0, 0.0]);
        assert_close(recompose(values, vectors), [[0.0; 3]; 3], 0.0);
    }

    #[test]
    fn a_known_decomposition_comes_back() {
        let matrix = [[4.0, 1.0, -2.0], [1.0, 2.0, 0.0], [-2.0, 0.0, 3.0]];

        let (values, vectors) = symmetric_eigen(matrix);

        assert!(values[0] >= values[1] && values[1] >= values[2]);
        assert_close(recompose(values, vectors), matrix, 1e-12);
    }

    #[test]
    fn the_eigenvectors_are_orthonormal() {
        let (_, v) = symmetric_eigen([[4.0, 1.0, -2.0], [1.0, 2.0, 0.0], [-2.0, 0.0, 3.0]]);

        for i in 0..3 {
            for j in 0..3 {
                let dot: f64 = (0..3).map(|r| v[r][i] * v[r][j]).sum();
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!((dot - expected).abs() < 1e-12, "columns {i},{j}: {dot}");
            }
        }
    }

    #[test]
    fn each_eigenpair_satisfies_its_own_equation() {
        let matrix = [[4.0, 1.0, -2.0], [1.0, 2.0, 0.0], [-2.0, 0.0, 3.0]];
        let (values, v) = symmetric_eigen(matrix);

        for k in 0..3 {
            for i in 0..3 {
                let av: f64 = (0..3).map(|j| matrix[i][j] * v[j][k]).sum();
                assert!(
                    (av - values[k] * v[i][k]).abs() < 1e-12,
                    "eigenpair {k}, component {i}: {} against {}",
                    av,
                    values[k] * v[i][k]
                );
            }
        }
    }

    #[test]
    fn two_equal_eigenvalues_still_decompose() {
        // The case a closed-form solution through the characteristic cubic
        // handles worst, and the reason this is Jacobi.
        let matrix = recompose(
            [5.0, 2.0, 2.0],
            [[0.0, 1.0, 0.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]],
        );

        let (values, vectors) = symmetric_eigen(matrix);

        assert!((values[0] - 5.0).abs() < 1e-12);
        assert!((values[1] - 2.0).abs() < 1e-12);
        assert!((values[2] - 2.0).abs() < 1e-12);
        assert_close(recompose(values, vectors), matrix, 1e-12);
    }

    #[test]
    fn a_scatter_matrix_of_coplanar_points_has_a_zero_eigenvalue() {
        // Points in the z = 0 plane: no scatter along z, so the smallest
        // eigenvalue vanishes and its eigenvector is the plane's normal. This
        // is the whole use this decomposition is put to.
        let points = [[1.0, 0.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, -2.0, 0.0]];

        let mut scatter = [[0.0; 3]; 3];
        for p in &points {
            for i in 0..3 {
                for j in 0..3 {
                    scatter[i][j] += p[i] * p[j];
                }
            }
        }

        let (values, vectors) = symmetric_eigen(scatter);

        assert!(values[2].abs() < 1e-12, "smallest eigenvalue {}", values[2]);
        assert!((vectors[2][2].abs() - 1.0).abs() < 1e-12, "normal is not along z");
    }

    #[test]
    fn a_very_flat_matrix_keeps_its_ordering() {
        // Fifteen orders of magnitude between the largest and smallest, which
        // is the range a nearly planar point cloud produces.
        let matrix = recompose(
            [1.0e6, 1.0e3, 1.0e-9],
            [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        );

        let (values, _) = symmetric_eigen(matrix);

        assert!(values[0] > values[1] && values[1] > values[2]);
        assert!((values[0] - 1.0e6).abs() < 1.0e-3);
        assert!(values[2] >= 0.0, "a scatter matrix has no negative eigenvalue");
    }
}
