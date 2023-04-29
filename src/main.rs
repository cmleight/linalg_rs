use num::{Zero, Float, One};

use std::fmt;
use std::ops::{Add, Div, Mul, Sub};
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
    } else {m};
    return res as usize
}

impl<T: IntFloat> Matrix<T> {
    pub fn laplace_determinant(&self) -> T {
        return (0..self.rows).map(|i| {
            (0..self.cols).map(|j| {
                (self.data[modulo_euclidean((i + j) as i64, self.rows as i64) + j * self.cols], if i + 1 < self.rows {
                    self.data[modulo_euclidean(i as i64 + 1 - j as i64, self.rows as i64) + j * self.cols]
                } else {
                    T::zero()
                })
            }).reduce(|(add_0, sub_0), (add_1, sub_1)| (add_0  * add_1, sub_0 * sub_1))
        })
        .map(|term| match term {
            Some((add, sub)) => add - sub,
            None => T::zero(),
        })
        .reduce(|a, b| a + b)
        .unwrap()
    }
}

fn main() {
    let a = Matrix {
        rows: 3,
        cols: 3,
        data: vec![1., 2., 2.,
                   2., 0., -1.,
                   -2., 1., 3.],
    };
    println!("determinant: {}", a.laplace_determinant())

}
