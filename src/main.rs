use num::{Float, One, Zero};

use std::fmt;
use std::ops::{Add, Div, Mul, Sub, Index, IndexMut};
use std::ops::{AddAssign, DivAssign, MulAssign, RemAssign, SubAssign};

/// use crate::ScalarOperand;

/// Elements that support linear algebra operations.
///
/// `'static` for type-based specialization, `Copy` so that they don't need move
/// semantics or destructors, and the rest are numerical traits.
pub trait LinalgScalar:
    'static
    + Copy
    + Zero
    + One
    + Add<Output = Self>
    + Sub<Output = Self>
    + Mul<Output = Self>
    + Div<Output = Self>
{
}

impl<T> LinalgScalar for T where
    T: 'static
        + Copy
        + Zero
        + One
        + Add<Output = T>
        + Sub<Output = T>
        + Mul<Output = T>
        + Div<Output = T>
{
}

impl IntFloat for f32 {}
impl IntFloat for f64 {}

pub trait IntFloat:
    Float
    + AddAssign
    + SubAssign
    + MulAssign
    + DivAssign
    + RemAssign
    + fmt::Display
    + fmt::Debug
    + fmt::LowerExp
    + fmt::UpperExp
    + LinalgScalar
    + Send
    + Sync
{
}

#[derive(Debug, PartialEq, Eq)]
struct Matrix<T: IntFloat> {
    rows: usize,
    cols: usize,
    data: Vec<T>,
}

impl<T: IntFloat> fmt::Display for Matrix<T> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let mut display_vec = String::with_capacity(self.rows * self.cols * 3);
        display_vec.push('[');
        display_vec.push('\n');
        (0..self.cols).for_each(|j| {
            display_vec.push('\t');
            display_vec.push('\t');
            display_vec.push('[');
            (0..self.rows)
                .for_each(|i| {
                    display_vec.push_str(&self[&j][i].to_string());
                    display_vec.push(',');
                });
            display_vec.pop();
            display_vec.push(']');
            display_vec.push(',');
            display_vec.push('\n');
        });
        display_vec.pop();
        display_vec.push('\n');
        display_vec.push('\t');
        display_vec.push(']');
        return write!(
            f,
            "{{\n\trows: {},\n\tcols: {},\n\tdata: {}\n}}",
            self.rows, self.cols, display_vec
        );
    }
}

impl<T:IntFloat> Index<&'_ usize> for Matrix<T> {
    type Output = [T];
    fn index(&self, i: &usize) -> &[T] {
        return &self.data[(i*self.rows)..(i*self.rows+self.cols)];
    }
}

impl<T:IntFloat> IndexMut<&'_ usize> for Matrix<T> {
    fn index_mut(&mut self, i: &'_ usize) -> &mut Self::Output {
        return &mut self.data[(i*self.rows)..(i*self.rows+self.cols)];
    }
}

impl<T: IntFloat> Matrix<T> {
    fn swap_rows(&mut self, i: usize, j: usize) {
        for k in 0..self.cols {
            self.data.swap(i * self.rows + k, j * self.rows + k)
        }
    }

    fn deep_copy(&self) -> Matrix<T> {
        return Matrix{
            rows: self.rows,
            cols: self.cols,
            // This deep copies the vector
            data: self.data.to_vec(),
        }
    }

    /// Below works for matricies <= rank 3
    /// use gaussian_elimination for rank >3
    pub fn laplace_determinant(&self) -> T {
        return (0..self.rows)
            .map(|i| {
                    (0..self.cols).map(|j| (
                            self[&j][(j + i) % self.cols],
                            self[&j][(self.cols + i - j) % self.cols]
                        )).reduce(|(a1, a2), (b1, b2)| {
                        (a1 * b1, a2 * b2)
                    }).unwrap()
            })
            .map(|(a, b)| a - b)
            .reduce(|a, b| a + b)
            .unwrap();
    }

    pub fn eigenvectors(&self) -> Option<Matrix<T>> {
        let a = self.deep_copy();
        if self.rows != self.cols { return None; }
        (0..self.rows);
        return Some(a);
    }

    pub fn eigenvalues(&self) -> Vec<T> {
        todo!()
    }


    /// Calculates the upper and lower triangular matricies from the source
    /// return the upper then lower matrix
    pub fn gaussian_elimination(&self) -> (bool, Matrix<T>) {
        let rows = self.rows;
        let cols = self.cols;
        let mut a = self.deep_copy();
        let mut det_sign_pos = true;
        let mut norm_col = 0;
        let mut norm_row = 0;

        while norm_col < cols && norm_row < rows {
            if a[&norm_col][norm_row].is_zero() {
                let i_max = (norm_col..rows)
                    .map(|i| (i, a[&i][norm_row]))
                    .reduce(|(x, a), (y, b)| if a.abs() >= b.abs() {(x, a)} else {(y, b)})
                    .map(|(x, _)| x)
                    .unwrap_or(norm_row);
                if a[&norm_col][i_max].is_zero() {
                    norm_col += 1;
                    continue;
                }
                if i_max != norm_row {
                    a.swap_rows(i_max, norm_row);
                    det_sign_pos = !det_sign_pos;
                }

            }
            ((norm_row + 1)..rows).for_each(|i| {
                let factor = a[&i][norm_col] / a[&norm_row][norm_col];
                a[&i][norm_col] = T::zero();
                ((norm_col + 1)..cols).for_each(|j| {
                    a[&i][j] = a[&i][j] - a[&norm_row][j] * factor;
                });
            });
            norm_col += 1;
            norm_row += 1;
        };

        return (det_sign_pos, a);
    }

    pub fn gaussian_determinant(&self) -> (Matrix<T>, T) {
        let (det_sign_pos, a) = self.gaussian_elimination();
        let base_det = (0..self.cols)
            .map(|i| a[&i][i])
            .reduce(|x, y| x * y)
            .unwrap();
        if !det_sign_pos {
            return (a, -base_det);
        } else {
            return (a, base_det);
        }
    }
}

fn main() { }

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_nonzero_laplace_determinant() {
        assert_eq!(Matrix {
            rows: 3,
            cols: 3,
            data: vec![1., 2., 2., 2., 0., -1., -2., 1., 3.],
        }.laplace_determinant(), -3.)
    }

    #[test]
    fn test_zero_laplace_determinant() {
        assert_eq!(Matrix {
            rows: 3,
            cols: 3,
            data: vec![1., 2., 3., 4., 5., 6., 7., 8., 9.],
        }.laplace_determinant(), 0.)
    }

    #[test]
    fn test_nonzero_gausian_elimination_determinant() {
        let res = Matrix {
            rows: 4,
            cols: 4,
            data: vec![1.,3.,5.,9.,1.,3.,1.,7.,4.,3.,9.,7.,5.,2.,0.,9.],
        }.gaussian_determinant();
        assert_eq!(res.1, -376.);
        assert_eq!(res.0, Matrix{
            rows: 4,
            cols: 4,
            data: vec![1.,3.,5.,9.,0.,-13.,-25.,-36.,0.,0.,6.307692307692307,-4.076923076923077,0.,0.,0.,-4.585365853658537],
        });
    }
}
