use std::cmp::{PartialEq, min};
use std::fmt::Display;
use std::f64::consts::PI;
use std::ops::{Add, Sub, Mul, Div};

/// Represents complex numbers with real and imaginary components.
/// In addition, an extra parameter--precision--is included for
/// displaying purposes.
#[derive(Debug, PartialEq, Clone, Copy)]
pub struct Complex {
    /// Real component
    pub real: f64,
    /// Imaginary component
    pub imag: f64,
    /// Number precision. Used for displaying purposes.
    pub prec: usize,
}

impl Complex {
    /// Construct a new complex number with a default precision of 3.
    /// 
    /// # Arguments
    /// * `real` - The real component
    /// * `imag` - The imaginary component
    /// 
    /// # Examples
    /// ```
    /// use sdomain_test::complex::Complex;
    /// 
    /// assert_eq!("1.000e0 + -2.000e0i", format!("{}", Complex::new(1.0, -2.0)));
    /// ```
    pub fn new(real: f64, imag: f64) -> Self {
        Self{real, imag, prec:3}
    }

    /// Construct a new complex number and set the precision.
    /// 
    /// # Arguments
    /// * `real` - The real component
    /// * `imag` - The imaginary component
    /// * `prec` - Precision
    /// 
    /// # Examples
    /// ```
    /// use sdomain_test::complex::Complex;
    /// 
    /// assert_eq!("1.0e0 + -2.0e0i", format!("{}", Complex::new_prec(1.0, -2.0, 1)));
    /// ```
    pub fn new_prec(real: f64, imag: f64, prec: usize) -> Self{
        Self{real, imag, prec}
    }

    /// Calculate the magnitude
    ///
    /// # Examples
    /// ```
    /// use sdomain_test::complex::Complex;
    /// 
    /// let mag = Complex::new(3.0, 4.0).mag();
    /// assert_eq!("5.0", format!("{:.1}", mag));
    /// ```
    pub fn mag(&self) -> f64 {
        f64::sqrt(f64::powf(self.real, 2.0) + f64::powf(self.imag, 2.0))
    }

    /// Calculate the magnitude and convert to decibels
    ///
    /// # Examples
    /// ```
    /// use sdomain_test::complex::Complex;
    /// 
    /// let mut mag_20log = Complex::new(3.0, 4.0).mag_20log();
    /// assert_eq!("13.979", format!("{:.3}", mag_20log));
    /// 
    /// mag_20log = Complex::new(1.0, 0.0).mag_20log();
    /// assert_eq!("0.000", format!("{:.3}", mag_20log));
    /// ```
    pub fn mag_20log(&self) -> f64 {
        20.0*self.mag().log10()
    }

    /// Calculate the phase in terms of radians
    /// 
    /// # Examples
    /// ```
    /// use sdomain_test::complex::Complex;
    /// 
    /// let phase = Complex::new(3.0, 4.0).phase();
    /// assert_eq!("0.927", format!("{:.3}", phase));
    /// ```
    pub fn phase(&self) -> f64 {
        f64::atan(self.imag / self.real)
    }

    /// Calculate the phase in terms of degrees
    /// 
    /// # Examples
    /// ```
    /// use sdomain_test::complex::Complex;
    /// 
    /// let phase_deg = Complex::new(3.0, 4.0).phase_deg();
    /// assert_eq!("53.130", format!("{:.3}", phase_deg));
    /// ```
    pub fn phase_deg(&self) -> f64 {
        self.phase()*180.0/PI
    }

    /// Get the complex conjugate
    /// 
    /// # Examples
    /// ```
    /// use sdomain_test::complex::Complex;
    /// 
    /// let c = Complex::new_prec(3.0, 4.0, 1);
    /// assert_eq!("3.0e0 + -4.0e0i", format!("{:}", c.conjugate()));
    /// ```
    pub fn conjugate(&self) -> Self {
        Self{real:self.real, imag:-self.imag, prec:self.prec}
    }
}

impl Display for Complex {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{0:.2$e} + {1:.2$e}i", self.real, self.imag, self.prec)
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
        assert_eq!("4.5000e0 + 6.8000e0i", format!("{}", Complex{real:4.500, imag:6.8, prec:4}));
    }

    #[test]
    fn conjugate_test() {
        assert_eq!("4.5000e0 + -6.8000e0i", format!("{}", Complex{real:4.500, imag:6.8, prec:4}.conjugate()));
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
        assert_eq!("5.00e0 + 1.10e1i", format!("{sum}"))
    }

    #[test]
    fn sub() {
        let (lhs, rhs) = gen_vals();
        let difference = lhs - rhs;
        assert_eq!("1.00e0 + -3.00e0i", format!("{difference}"))
    }

    #[test]
    fn mul() {
        let (lhs, rhs) = gen_vals();
        let product = lhs * rhs;
        assert_eq!("-2.20e1 + 2.90e1i", format!("{product}"))
    }

    #[test]
    fn div() {
        let (lhs, rhs) = gen_vals();
        let quotient = lhs / rhs;
        assert_eq!("6.42e-1 + -2.45e-1i", format!("{quotient}"))
    }
}