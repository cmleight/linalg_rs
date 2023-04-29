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

impl<T: IntFloat> Matrix<T> {
    pub fn determinant(&self) {
        (0..self.rows).map(|i| {
            
        })
    }
}

fn main() {
    let a = Matrix {
        rows: 3,
        cols: 3,
        data: vec![1., 2., 3.,
                   4., 5., 6.,
                   7., 8., 9.],
    };

}
