use std::cmp::{PartialEq, min};
use std::fmt::Display;
use std::f64::consts::PI;
use std::ops::{Add, Sub, Mul, Div};

#[derive(Debug, PartialEq, Clone, Copy)]
pub struct Complex {
    pub real: f64,
    pub imag: f64,
    pub prec: usize,
}

impl Complex {
    pub fn new(real: f64, imag: f64) -> Self {
        Self{real, imag, prec:3}
    }

    pub fn new_prec(real: f64, imag: f64, prec: usize) -> Self{
        Self{real, imag, prec}
    }

    pub fn mag(&self) -> f64 {
        f64::sqrt(f64::powf(self.real, 2.0) + f64::powf(self.imag, 2.0))
    }

    pub fn mag_20log(&self) -> f64 {
        20.0*self.mag().log10()
    }

    pub fn phase(&self) -> f64 {
        f64::atan(self.imag / self.real)
    }

    pub fn phase_deg(&self) -> f64 {
        self.phase()*180.0/PI
    }

    pub fn conjugate(&self) -> Self {
        Self{real:self.real, imag:-self.imag, prec:self.prec}
    }
}

impl Display for Complex {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{0:.2$} + {1:.2$}i", self.real, self.imag, self.prec)
    }
}

impl Add for Complex {
    type Output = Complex;

    fn add(self, rhs: Self) -> Self {
        Self{
            real: self.real + rhs.real,
            imag: self.imag + rhs.imag,
            prec: min(self.prec, rhs.prec)
        }
    }
}

impl Sub for Complex {
    type Output = Complex;

    fn sub(self, rhs: Self) -> Self {
        Self{
            real: self.real - rhs.real,
            imag: self.imag - rhs.imag,
            prec: min(self.prec, rhs.prec)
        }
    }
}

impl Mul for Complex {
    type Output = Complex;

    fn mul(self, rhs: Self) -> Self {
        Self{
            real: self.real*rhs.real - self.imag*rhs.imag,
            imag: self.real*rhs.imag + self.imag*rhs.real,
            prec: min(self.prec, rhs.prec)
        }
    }
}

impl Div for Complex {
    type Output = Complex;

    fn div(self, rhs: Self) -> Self {
        let divisor = rhs.real*rhs.real + rhs.imag*rhs.imag;
        Self{
            real: (self.real*rhs.real + self.imag*rhs.imag) / divisor,
            imag: (-self.real*rhs.imag + self.imag*rhs.real) / divisor,
            prec: min(self.prec, rhs.prec)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn new_test() {
        assert_eq!(Complex::new(4.0, 5.0), Complex{real:4.0, imag:5.0, prec:3});
    }

    #[test]
    fn new_prec_test() {
        assert_eq!(Complex::new_prec(4.0, 5.0, 67), Complex{real:4.0, imag:5.0, prec:67});
    }

    #[test]
    fn mag_test() {
        assert_eq!(Complex::new(4.0, 3.0).mag(), 5.0);
    }

    #[test]
    fn phase_test() {
        assert_eq!(Complex::new(1.0, 0.0).phase(), 0.0);
    }

    #[test]
    fn phase_gt1_test() {
        assert_eq!(Complex::new(2.0, 0.0).phase(), 0.0);
    }

    #[test]
    fn phase_div0_test() {
        assert_eq!(Complex::new(0.0, 1.0).phase(), PI/2.0);
    }

    #[test]
    fn phase_deg_test() {
        assert_eq!(Complex::new(0.0, 1.0).phase_deg(), 90.0);
    }

    #[test]
    fn fmt_test() {
        assert_eq!("4.5000 + 6.8000i", format!("{}", Complex{real:4.500, imag:6.8, prec:4}));
    }

    #[test]
    fn conjugate_test() {
        assert_eq!("4.5000 + -6.8000i", format!("{}", Complex{real:4.500, imag:6.8, prec:4}.conjugate()));
    }
}

#[cfg(test)]
mod arith_tests {
    use super::*;

    fn gen_vals() -> (Complex, Complex) {
        (
            Complex::new(3.0, 4.0),
            Complex::new_prec(2.0, 7.0, 2)
        )
    }

    #[test]
    fn add() {
        let (lhs, rhs) = gen_vals();
        let sum = lhs + rhs;
        assert_eq!("5.00 + 11.00i", format!("{sum}"))
    }

    #[test]
    fn sub() {
        let (lhs, rhs) = gen_vals();
        let difference = lhs - rhs;
        assert_eq!("1.00 + -3.00i", format!("{difference}"))
    }

    #[test]
    fn mul() {
        let (lhs, rhs) = gen_vals();
        let product = lhs * rhs;
        assert_eq!("-22.00 + 29.00i", format!("{product}"))
    }

    #[test]
    fn div() {
        let (lhs, rhs) = gen_vals();
        let quotient = lhs / rhs;
        assert_eq!("0.64 + -0.25i", format!("{quotient}"))
    }
}