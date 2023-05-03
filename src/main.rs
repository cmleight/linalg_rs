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
///    + ScalarOperand
#[derive(Debug)]
struct Matrix<T: IntFloat> {
    rows: usize,
    cols: usize,
    data: Vec<T>,
}

fn modulo_euclidean(a: i64, b: i64) -> usize {
    let m = a % b;
    let res = if m < 0 {
        if b < 0 {
            m - b
        } else {
            m + b
        }
    } else {
        m
    };
    return res as usize;
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
                    display_vec.push_str(&self.data[j + i * self.rows].to_string());
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
    fn deep_copy(&self) -> Matrix<T> {
        return Matrix{
            rows: self.rows,
            cols: self.cols,
            // This deep copies the vector
            data: self.data.to_vec(),
        }
    }

    pub fn laplace_determinant(&self) -> T {
        return (0..self.rows)
            .map(|i| {
                    (0..self.cols).map(|j| (
                        self[&j][(j + i) % self.cols],
                        self[&j][(self.cols + i - j) % self.cols]
                    )).reduce(|(a1, a2), (b1, b2)| (a1 * a2, b1 * b2)).unwrap()
            })
            .map(|(a, b)| a - b)
            .reduce(|a, b| a + b)
            .unwrap();
    }

    /// Calculates the upper and lower triangular matricies from the source
    /// return the upper then lower matrix
    pub fn gaussian_elimination(&self) -> Matrix<T> {
        let rows = self.rows;
        let cols = self.cols;
        let mut a = self.deep_copy();

        (0..rows).for_each(|j|{
            // pivot is not zero
            // TODO: swap rows when zero
            if a[&j][j] != T::zero() {
                (j..cols).for_each(|i|{

                })
            }
        });

        return a
    }

    //    pub fn lu_determinants(&self) -> T {
    //    }
}

fn main() {
    let a = Matrix {
        rows: 3,
        cols: 3,
        data: vec![1., 2., 2., 2., 0., -1., -2., 1., 3.],
    };
    println!("base_matrix: {}", a);
    println!("determinant: {}", a.laplace_determinant());
}
