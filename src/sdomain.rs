use std::{
    fmt::Display,
    f64::consts::PI,
    ops::{
        Neg,
        Add, Sub, Mul, Div,
        DivAssign, AddAssign, SubAssign, MulAssign
    },
    cmp::min,
    mem::swap, iter::Sum
};

use crate::complex::Complex;

/// Represents a model in the s-domain as a ratio of polynomials.
#[derive(Debug, PartialEq, Clone)]
pub struct Fs {
    /// Polynomial in terms of s in the numerator.
    pub numer: Vec<f64>,
    /// Polynomial in terms of s in the denominator.
    pub denom: Vec<f64>,
    /// Used for displaying values. Number of digits after the decimal in scientific notation.
    pub prec: usize,
}

impl Fs {
    /// Helper to calculate the output of a polynomial for a given frequency in radians.
    /// Returns a complex number.
    /// 
    /// # Arguments
    /// * `polynomial` - Collection holding the coefficients of the polynomial.
    /// * `rad` - Frequency in radians to calculate for.
    /// * `precision` - Precision to give to the resulting complex value.
    fn merge_rad(polynomial: &Vec<f64>, rad: f64, precision: usize) -> Complex {
        let mut real = 0.0;
        let mut imag = 0.0;
        let mut val = 1.0;
        for (idx, coeff) in polynomial.iter().enumerate() {
            if idx % 2 == 0 { real += coeff*val*(if idx % 4 == 0 {1.0} else {-1.0}); }
            else            { imag += coeff*val*(if (idx-1) % 4 == 0 {1.0} else {-1.0}); }
            val *= rad;
        }
        Complex { real, imag, prec: precision }
    }

    /// Helper to calculate the output of a polynomial for a given frequency in Hertz.
    /// Returns a complex number.
    /// 
    /// # Arguments
    /// * `polynomial` - Collection holding the coefficients of the polynomial.
    /// * `rad` - Frequency in Hertz to calculate for.
    /// * `precision` - Precision to give to the resulting complex value.
    fn merge_freq(polynomial: &Vec<f64>, freq: f64, precision: usize) -> Complex {
        Self::merge_rad(polynomial, freq*2.0*PI, precision)
    }

    /// Helper to generate superscripts for displaying exponents
    /// Returns a String object.
    /// 
    /// # Arguments
    /// * `exponent` - The exponent to convert to a superscript string.
    fn to_superscript(mut exponent: usize) -> String {
        if exponent == 0 {return "⁰".to_string();}
        let mut s = String::new();
        while exponent > 0 {
            let modulus = exponent % 10;
            exponent /= 10;
            s.insert(0, match modulus {
                0 => '⁰',
                1 => '¹',
                2 => '²',
                3 => '³',
                4 => '⁴',
                5 => '⁵',
                6 => '⁶',
                7 => '⁷',
                8 => '⁸',
                9 => '⁹',
                _ => panic!("modulus is unexpectedly outside of the range [0, 9]")
            });
        }
        s
    }

    /// Helper to format an s^n term for displaying.
    /// Returns a String object.
    /// 
    /// # Arguments
    /// * `exponent` - The power of s: s^exponent
    fn fmt_s_term(exponent: usize) -> String {
        match exponent {
            0 => String::from(""),
            1 => String::from("s"),
            _ => format!("s{}", Self::to_superscript(exponent)),
            // _ => format!("s^{exponent}"),
        }
    }

    /// Helper to format a polynomial for displaying.
    /// Returns a String object.
    /// 
    /// # Arguments
    /// * `polynomial` - Contains the coefficients of the polynomial.
    /// * `precision` - Determines how many digits to display after the decimal point.
    /// * `show_all` - If False, "s⁰" and "s¹" are shortened to "" and "s", respectively.
    fn to_str(polynomial: &Vec<f64>, precision: usize, show_all: bool) -> String {
        let terms: Vec<String> = polynomial.iter().enumerate().map(|(idx, val)| {
            if show_all {
                format!("{val:.0$e}s{1}", precision, Self::to_superscript(idx))
            } else {
                format!("{val:.0$e}{1}", precision, Self::fmt_s_term(idx))
            }
        }).collect();
        let mut toreturn = String::new();
        for term in terms {
            toreturn.push_str(&term);
            toreturn.push_str(" + ");
        }
        toreturn[..(toreturn.len()-3)].to_string()
    }

    /// Formats the numerator for displaying.
    /// Returns a String object.
    /// 
    /// # Arguments
    /// * `show_all` - If False, "s⁰" and "s¹" are shortened to "" and "s", respectively.
    /// 
    /// # Examples
    /// ```
    /// use sdomain_test::sdomain::Fs;
    /// 
    /// let fs = Fs{numer: vec![1.0, 2.0, 3.0], denom: vec![4.0, 5.0, 6.0], prec:1};
    /// assert_eq!(
    ///     "1.0e0 + 2.0e0s + 3.0e0s²".to_string(),
    ///     fs.to_str_numer(false)
    /// );
    /// assert_eq!(
    ///     "1.0e0s⁰ + 2.0e0s¹ + 3.0e0s²".to_string(),
    ///     fs.to_str_numer(true)
    /// );
    /// ```
    pub fn to_str_numer(&self, show_all: bool) -> String {
        Self::to_str(&self.numer, self.prec, show_all)
    }

    /// Formats the denominator for displaying.
    /// Returns a String object.
    /// 
    /// # Arguments
    /// * `show_all` - If False, "s⁰" and "s¹" are shortened to "" and "s", respectively.
    /// 
    /// # Examples
    /// ```
    /// use sdomain_test::sdomain::Fs;
    /// 
    /// let fs = Fs{numer: vec![1.0, 2.0, 3.0], denom: vec![4.0, 5.0, 6.0], prec:1};
    /// assert_eq!(
    ///     "4.0e0 + 5.0e0s + 6.0e0s²".to_string(),
    ///     fs.to_str_denom(false)
    /// );
    /// assert_eq!(
    ///     "4.0e0s⁰ + 5.0e0s¹ + 6.0e0s²".to_string(),
    ///     fs.to_str_denom(true)
    /// );
    /// ```
    pub fn to_str_denom(&self, show_all: bool) -> String {
        Self::to_str(&self.denom, self.prec, show_all)
    }

    /// Calculate the output of an s-domain model for a given frequency in radians.
    /// Returns a complex number.
    /// 
    /// # Arguments
    /// * `fs` - The s-domain model
    /// * `rad` - The frequency in radians to calculate for
    /// 
    /// # Examples
    /// ```
    /// use sdomain_test::sdomain::Fs;
    /// use sdomain_test::complex::Complex;
    /// 
    /// let fs = Fs{numer: vec![1.0, 2.0], denom: vec![-0.5, 2.5], prec: 1};
    /// assert_eq!("8.0e-1 + -5.6e-5i", format!("{}", fs.calculate_rad(10e3)));
    /// ```
    pub fn calculate_rad(&self, rad: f64) -> Complex {
        Fs::merge_rad(&self.numer, rad, self.prec) / Fs::merge_rad(&self.denom, rad, self.prec)
    }

    /// Calculate the output of an s-domain model for a given frequency in Hertz.
    /// Returns a complex number.
    /// 
    /// # Arguments
    /// * `fs` - The s-domain model
    /// * `freq` - The frequency in Hertz to calculate for
    /// 
    /// # Examples
    /// ```
    /// use sdomain_test::sdomain::Fs;
    /// use sdomain_test::complex::Complex;
    /// use std::f64::consts::PI;
    /// 
    /// let fs = Fs{numer: vec![1.0, 2.0], denom: vec![-0.5, 2.5], prec: 1};
    /// assert_eq!("8.0e-1 + -5.6e-5i", format!("{}", fs.calculate_freq(10e3/2.0/PI)));
    /// ```
    pub fn calculate_freq(&self, freq: f64) -> Complex {
        Fs::merge_freq(&self.numer, freq, self.prec) / Fs::merge_freq(&self.denom, freq, self.prec)
    }

    /// Multiply the model by some scaling value.
    /// Returns the same object for chaining function calls.
    /// 
    /// # Arguments
    /// * `scale` - The value to scale by.
    /// 
    /// # Examples
    /// ```
    /// use sdomain_test::sdomain::Fs;
    /// 
    /// let fs = Fs{numer: vec![1.0, 2.0], denom: vec![-0.5, 2.5], prec: 1};
    /// assert_eq!(
    ///     "( 3.0e0 + 6.0e0s ) / ( -5.0e-1 + 2.5e0s )",
    ///     format!("{}", fs.scale(3.0))
    /// );
    /// 
    /// let fs2 = Fs{numer: vec![1.0, 2.0], denom: vec![-0.5, 2.5], prec: 1};
    /// assert_eq!(
    ///     "( 1.0e0 + 2.0e0s ) / ( -5.0e-1 + 2.5e0s )",
    ///     format!("{}", fs2.scale(3.0).scale(1.0/3.0))
    /// );
    /// ```
    pub fn scale(mut self, scale: f64) -> Self {
        for term in self.numer.iter_mut() {
            *term *= scale;
        }
        self
    }

    /// The same as raising the model to a power of -1. Swaps the numerator and denominator.
    /// Returns the same object for chaining function calls.
    /// 
    /// # Examples
    /// ```
    /// use sdomain_test::sdomain::Fs;
    /// 
    /// let fs = Fs{numer: vec![1.0, 2.0], denom: vec![-0.5, 2.5], prec: 1};
    /// assert_eq!(
    ///     "( -5.0e-1 + 2.5e0s ) / ( 1.0e0 + 2.0e0s )",
    ///     format!("{}", fs.invert())
    /// );
    /// 
    /// let fs2 = Fs{numer: vec![1.0, 2.0], denom: vec![-0.5, 2.5], prec: 1};
    /// assert_eq!(
    ///     "( 1.0e0 + 2.0e0s ) / ( -5.0e-1 + 2.5e0s )",
    ///     format!("{}", fs2.invert().invert())
    /// );
    /// ```
    pub fn invert(mut self) -> Self {
        swap(&mut self.numer, &mut self.denom);
        self
    }

    /// Helper to calculate the sum of a collection of polynomials.
    /// Returns a polynomial representing the sum.
    /// 
    /// # Arguments
    /// * `bucket` - A collection of polynomials.
    fn accumulate<T>(bucket: T) -> Vec<f64> where
        T: IntoIterator<Item = Vec<f64>>,
        T::Item: IntoIterator,
        <T::Item as IntoIterator>::Item: Add,
    {
        let mut iter = bucket.into_iter();
        let mut sum = iter.next().unwrap();
        for mut factor in iter {
            if sum.len() > factor.len() {
                factor.resize(sum.len(), 0.0);
            } else {
                sum.resize(factor.len(), 0.0);
            }
            for idx in 0..sum.len() {
                sum[idx] += factor[idx];
            }
        }
        sum
    }

    /// Helper to handle raw multiplication between two polynomials.
    /// Returns a polynomial representing the product.
    /// 
    /// # Arguments
    /// * `lhs` - The first polynomial to multiply by
    /// * `rhs` - The second polynomial to multiply by
    fn multiply_helper(lhs: &Vec<f64>, rhs: &Vec<f64>) -> Vec<f64> {
        let mut bucket: Vec::<Vec::<f64>> = rhs.into_iter()
            .map(|rhs_coeff| {
                lhs.into_iter()
                    .map(|x| x*rhs_coeff)
                    .collect()
            })
            .collect();
        let mut resized_bucket = Vec::<Vec::<f64>>::with_capacity(bucket.len());
        for (idx, factor) in bucket.iter_mut().enumerate() {
            resized_bucket.push(factor.splice(0..0, vec![0.0; idx]).collect::<Vec::<f64>>());
        }
        Self::accumulate(bucket)
    }

    /// Helper to trim unneeded, leading zeroes from both the numerator and denominator.
    /// For example, (1 + 0s + 2s² + 0s³) becomes (1 + 0s + 2s²).
    fn trim_leading_zeroes(&mut self) {
        Self::trim_helper(&mut self.numer);
        Self::trim_helper(&mut self.denom);
    }

    /// Helper to trim unneeded, leading zeroes from a polynomial.
    /// For example, (1 + 0s + 2s² + 0s³) becomes (1 + 0s + 2s²).
    /// 
    /// # Arguments
    /// * `polynomial` - The polynomial to trim.
    fn trim_helper(polynomial: &mut Vec<f64>) {
        let mut leading_count = 0;
        for val in polynomial.iter().rev() {
            if format!("{val:.0e}").chars().nth(0).unwrap() == '0' {leading_count += 1;}
            else {break;}
        }
        // If the entire collection is 0.0, keep at least one
        leading_count = min(leading_count, polynomial.len()-1);
        polynomial.truncate(polynomial.len() - leading_count);
    }
}

impl Display for Fs {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let numer_str = self.to_str_numer(false);
        let denom_str = self.to_str_denom(false);
        write!(f, "( {0} ) / ( {1} )", numer_str, denom_str)
    }
}

impl Complex {
    /// Calculate the output of an s-domain model for a given frequency in radians.
    /// Returns a complex number.
    /// 
    /// # Arguments
    /// * `fs` - The s-domain model
    /// * `rad` - The frequency in radians to calculate for
    /// 
    /// # Examples
    /// ```
    /// use sdomain_test::sdomain::Fs;
    /// use sdomain_test::complex::Complex;
    /// 
    /// let fs = Fs{numer: vec![1.0, 2.0], denom: vec![-0.5, 2.5], prec: 1};
    /// assert_eq!("8.0e-1 + -5.6e-5i", format!("{}", Complex::from_rad(fs, 10e3)));
    /// ```
    pub fn from_rad(fs: Fs, rad: f64) -> Self {
        fs.calculate_rad(rad)
    }

    /// Calculate the output of an s-domain model for a given frequency in Hertz.
    /// Returns a complex number.
    /// 
    /// # Arguments
    /// * `fs` - The s-domain model
    /// * `freq` - The frequency in Hertz to calculate for
    /// 
    /// # Examples
    /// ```
    /// use sdomain_test::sdomain::Fs;
    /// use sdomain_test::complex::Complex;
    /// use std::f64::consts::PI;
    /// 
    /// let fs = Fs{numer: vec![1.0, 2.0], denom: vec![-0.5, 2.5], prec: 1};
    /// assert_eq!("8.0e-1 + -5.6e-5i", format!("{}", Complex::from_freq(fs, 10e3/2.0/PI)));
    /// ```
    pub fn from_freq(fs: Fs, freq: f64) -> Self {
        fs.calculate_freq(freq)
    }
}

impl Neg for Fs {
    type Output = Fs;

    fn neg(self) -> Self::Output {
        Fs{numer:self.numer.iter().map(|x| x*-1.0).collect(), denom:self.denom, prec:self.prec}
    }
}

impl<'b> Add<&'b Fs> for Fs {
    type Output = Fs;

    fn add(self, rhs: &'b Fs) -> Self::Output {
        let mut toreturn = self.clone();
        toreturn += rhs;
        toreturn
    }
}

impl<'b> Sub<&'b Fs> for Fs {
    type Output = Fs;

    fn sub(self, rhs: &'b Fs) -> Self::Output {
        self + &-rhs.clone()
    }
}

impl<'b> Mul<&'b Fs> for Fs {
    type Output = Fs;

    fn mul(self, rhs: &'b Fs) -> Self::Output {
        let mut toreturn = self.clone();
        toreturn *= rhs;
        toreturn
    }
}

impl<'b> Div<&'b Fs> for Fs {
    type Output = Fs;

    fn div(self, rhs: &'b Fs) -> Self::Output {
        let mut toreturn = self.clone();
        toreturn /= rhs;
        toreturn
    }
}

impl<'b> AddAssign<&'b Fs> for Fs {
    fn add_assign(&mut self, rhs: &'b Fs) {
        if self.denom == rhs.denom {
            self.numer = Fs::accumulate([self.numer.clone(), rhs.numer.clone()]);
            self.prec = min(self.prec, rhs.prec);
        } else {
            let lhs_denom = self.denom.clone();
            let rhs_denom = rhs.denom.clone();
            self.denom = Fs::multiply_helper(&self.denom, &rhs.denom);
            self.numer = Fs::accumulate([
                    Fs::multiply_helper(&self.numer, &rhs_denom),
                    Fs::multiply_helper(&rhs.numer, &lhs_denom)
                    ]);
            self.prec = min(self.prec, rhs.prec);
        }
        self.trim_leading_zeroes();
    }
}

impl<'b> SubAssign<&'b Fs> for Fs {
    fn sub_assign(&mut self, rhs: &'b Fs) {
        *self += &-rhs.clone();
    }
}

impl<'b> MulAssign<&'b Fs> for Fs {
    fn mul_assign(&mut self, rhs: &'b Fs) {
        self.numer = Fs::multiply_helper(&self.numer, &rhs.numer);
        self.denom = Fs::multiply_helper(&self.denom, &rhs.denom);
        self.prec = min(self.prec, rhs.prec);
    }
}

impl<'b> DivAssign<&'b Fs> for Fs {
    fn div_assign(&mut self, rhs: &'b Fs) {
        self.numer = Fs::multiply_helper(&self.numer, &rhs.denom);
        self.denom = Fs::multiply_helper(&self.denom, &rhs.numer);
        self.prec = min(self.prec, rhs.prec);
    }
}

impl Sum for Fs {
    fn sum<I: Iterator<Item = Self>>(mut iter: I) -> Self {
        match iter.next() {
            Some(acc) => iter.fold(acc, |a, fs| a + &fs),
            None => gen::zero_prec(12)
        }
    }
}

pub mod gen {
    use super::Fs;

    /// Generate a model in the s-domain equivalent to 0.0.
    /// 
    /// # Arguments
    /// * `prec` - Used for displaying values. Number of digits after the decimal in scientific notation.
    /// 
    /// # Examples
    /// ```
    /// use sdomain_test::sdomain::gen;
    /// 
    /// assert_eq!("( 0.000e0 ) / ( 1.000e0 )", format!("{}", gen::zero_prec(3)));
    /// ```
    pub fn zero_prec(prec: usize) -> Fs {Fs{numer:vec![0.0], denom:vec![1.0], prec}}
    /// Generate a model in the s-domain equivalent to 1.0.
    /// 
    /// # Arguments
    /// * `prec` - Used for displaying values. Number of digits after the decimal in scientific notation.
    /// 
    /// # Examples
    /// ```
    /// use sdomain_test::sdomain::gen;
    /// 
    /// assert_eq!("( 1.000e0 ) / ( 1.000e0 )", format!("{}", gen::unit_prec(3)));
    /// ```
    pub fn unit_prec(prec: usize) -> Fs {Fs{numer:vec![1.0], denom:vec![1.0], prec}}
    /// Generate a model in the s-domain equivalent to s.
    /// 
    /// # Arguments
    /// * `prec` - Used for displaying values. Number of digits after the decimal in scientific notation.
    /// 
    /// # Examples
    /// ```
    /// use sdomain_test::sdomain::gen;
    /// 
    /// assert_eq!("( 0.000e0 + 1.000e0s ) / ( 1.000e0 )", format!("{}", gen::s_prec(3)));
    /// ```
    pub fn s_prec(prec: usize) -> Fs {Fs{numer:vec![0.0, 1.0], denom:vec![1.0], prec}}
    /// Generate a model in the s-domain equivalent to 1/s.
    /// 
    /// # Arguments
    /// * `prec` - Used for displaying values. Number of digits after the decimal in scientific notation.
    /// 
    /// # Examples
    /// ```
    /// use sdomain_test::sdomain::gen;
    /// 
    /// assert_eq!("( 1.000e0 ) / ( 0.000e0 + 1.000e0s )", format!("{}", gen::step_prec(3)));
    /// ```
    pub fn step_prec(prec: usize) -> Fs {Fs{numer:vec![1.0], denom:vec![0.0, 1.0], prec}}
    /// Generate an impedance model in the s-domain representing an ideal resistor.
    /// 
    /// # Arguments
    /// * `r` - Resistance in Ohms.
    /// * `prec` - Used for displaying values. Number of digits after the decimal in scientific notation.
    /// 
    /// # Examples
    /// ```
    /// use sdomain_test::sdomain::gen;
    /// 
    /// assert_eq!("( 2.000e1 ) / ( 1.000e0 )", format!("{}", gen::resistor_prec(20.0, 3)));
    /// ```
    pub fn resistor_prec(r: f64, prec: usize) -> Fs {Fs{numer:vec![r], denom:vec![1.0], prec}}
    /// Generate an impedance model in the s-domain representing an ideal capacitor.
    /// 
    /// # Arguments
    /// * `c` - Capacitance in Farads.
    /// * `prec` - Used for displaying values. Number of digits after the decimal in scientific notation.
    /// 
    /// # Examples
    /// ```
    /// use sdomain_test::sdomain::gen;
    /// 
    /// assert_eq!("( 1.000e0 ) / ( 0.000e0 + 3.300e-6s )", format!("{}", gen::capacitor_prec(3.3e-6, 3)));
    /// ```
    pub fn capacitor_prec(c: f64, prec: usize) -> Fs {Fs{numer:vec![1.0], denom:vec![0.0, c], prec}}
    /// Generate an impedance model in the s-domain representing an ideal inductor.
    /// 
    /// # Arguments
    /// * `l` - Inductance in Henries.
    /// * `prec` - Used for displaying values. Number of digits after the decimal in scientific notation.
    /// 
    /// # Examples
    /// ```
    /// use sdomain_test::sdomain::gen;
    /// 
    /// assert_eq!("( 0.000e0 + 3.300e-6s ) / ( 1.000e0 )", format!("{}", gen::inductor_prec(3.3e-6, 3)));
    /// ```
    pub fn inductor_prec(l: f64, prec: usize) -> Fs {Fs{numer:vec![0.0, l], denom:vec![1.0], prec}}

    /// Generate an impedance model in the s-domain representing an RL series circuit or an inductor with DCR.
    /// 
    /// # Arguments
    /// * `r` - Resistance in Ohms.
    /// * `l` - Inductance in Henries.
    /// * `prec` - Used for displaying values. Number of digits after the decimal in scientific notation.
    /// 
    /// # Examples
    /// ```
    /// use sdomain_test::sdomain::gen;
    /// 
    /// assert_eq!("( 4.000e-2 + 3.300e-6s ) / ( 1.000e0 )", format!("{}", gen::rl_prec(0.04, 3.3e-6, 3)));
    /// ```
    pub fn rl_prec(r: f64, l: f64, prec: usize) -> Fs {Fs{numer:vec![r, l], denom:vec![1.0], prec}}
    /// Generate an impedance model in the s-domain representing an RCL series circuit or a capacitor with ESR and ESL.
    /// 
    /// # Arguments
    /// * `r` - Resistance in Ohms.
    /// * `c` - Capacitance in Farads.
    /// * `l` - Inductance in Henries.
    /// * `prec` - Used for displaying values. Number of digits after the decimal in scientific notation.
    /// 
    /// # Examples
    /// ```
    /// use sdomain_test::sdomain::gen;
    /// 
    /// assert_eq!("( 1.000e0 + 6.600e-9s + 4.950e-15s² ) / ( 0.000e0 + 3.300e-6s )", format!("{}", gen::rcl_prec(2e-3, 3.3e-6, 1.5e-9, 3)));
    /// ```
    pub fn rcl_prec(r: f64, c: f64, l: f64, prec: usize) -> Fs {Fs{numer:vec![1.0, r*c, l*c], denom:vec![0.0, c], prec}}



    /// Generate a model in the s-domain equivalent to 0.0.
    /// 
    /// # Examples
    /// ```
    /// use sdomain_test::sdomain::gen;
    /// 
    /// assert_eq!("( 0.000e0 ) / ( 1.000e0 )", format!("{}", gen::zero()));
    /// ```
    pub fn zero() -> Fs {zero_prec(3)}
    /// Generate a model in the s-domain equivalent to 1.0.
    /// 
    /// # Examples
    /// ```
    /// use sdomain_test::sdomain::gen;
    /// 
    /// assert_eq!("( 1.000e0 ) / ( 1.000e0 )", format!("{}", gen::unit()));
    /// ```
    pub fn unit() -> Fs {unit_prec(3)}
    /// Generate a model in the s-domain equivalent to s.
    /// 
    /// # Examples
    /// ```
    /// use sdomain_test::sdomain::gen;
    /// 
    /// assert_eq!("( 0.000e0 + 1.000e0s ) / ( 1.000e0 )", format!("{}", gen::s()));
    /// ```
    pub fn s() -> Fs {s_prec(3)}
    /// Generate a model in the s-domain equivalent to 1/s.
    /// 
    /// # Examples
    /// ```
    /// use sdomain_test::sdomain::gen;
    /// 
    /// assert_eq!("( 1.000e0 ) / ( 0.000e0 + 1.000e0s )", format!("{}", gen::step()));
    /// ```
    pub fn step() -> Fs {step_prec(3)}
    /// Generate an impedance model in the s-domain representing an ideal resistor.
    /// 
    /// # Arguments
    /// * `r` - Resistance in Ohms.
    /// 
    /// # Examples
    /// ```
    /// use sdomain_test::sdomain::gen;
    /// 
    /// assert_eq!("( 2.000e1 ) / ( 1.000e0 )", format!("{}", gen::resistor(20.0)));
    /// ```
    pub fn resistor(r: f64) -> Fs {resistor_prec(r, 3)}
    /// Generate an impedance model in the s-domain representing an ideal capacitor.
    /// 
    /// # Arguments
    /// * `c` - Capacitance in Farads.
    /// 
    /// # Examples
    /// ```
    /// use sdomain_test::sdomain::gen;
    /// 
    /// assert_eq!("( 1.000e0 ) / ( 0.000e0 + 3.300e-6s )", format!("{}", gen::capacitor(3.3e-6)));
    /// ```
    pub fn capacitor(c: f64) -> Fs {capacitor_prec(c, 3)}
    /// Generate an impedance model in the s-domain representing an ideal inductor.
    /// 
    /// # Arguments
    /// * `l` - Inductance in Henries.
    /// 
    /// # Examples
    /// ```
    /// use sdomain_test::sdomain::gen;
    /// 
    /// assert_eq!("( 0.000e0 + 3.300e-6s ) / ( 1.000e0 )", format!("{}", gen::inductor(3.3e-6)));
    /// ```
    pub fn inductor(l: f64) -> Fs {inductor_prec(l, 3)}

    /// Generate an impedance model in the s-domain representing an RL series circuit or an inductor with DCR.
    /// 
    /// # Arguments
    /// * `r` - Resistance in Ohms.
    /// * `l` - Inductance in Henries.
    /// 
    /// # Examples
    /// ```
    /// use sdomain_test::sdomain::gen;
    /// 
    /// assert_eq!("( 4.000e-2 + 3.300e-6s ) / ( 1.000e0 )", format!("{}", gen::rl(0.04, 3.3e-6)));
    /// ```
    pub fn rl(r: f64, l: f64) -> Fs {rl_prec(r, l, 3)}
    /// Generate an impedance model in the s-domain representing an RCL series circuit or a capacitor with ESR and ESL.
    /// 
    /// # Arguments
    /// * `r` - Resistance in Ohms.
    /// * `c` - Capacitance in Farads.
    /// * `l` - Inductance in Henries.
    /// 
    /// # Examples
    /// ```
    /// use sdomain_test::sdomain::gen;
    /// 
    /// assert_eq!("( 1.000e0 + 6.600e-9s + 4.950e-15s² ) / ( 0.000e0 + 3.300e-6s )", format!("{}", gen::rcl(2e-3, 3.3e-6, 1.5e-9)));
    /// ```
    pub fn rcl(r: f64, c: f64, l: f64) -> Fs {rcl_prec(r, c, l, 3)}
}

/// Combine two models assuming they are in parallel.
/// Returns a combined model.
/// ```text
/// 
///           ┌───┴───┐
///        ┌──┴─┐  ┌──┴─┐
/// H(s) = │F(s)│  │G(s)│
///        └──┬─┘  └──┬─┘
///           └───┬───┘
/// 
/// ```
/// 
/// # Arguments
/// * `fs1` - One of the models in parallel
/// * `fs2` - One of the models in parallel
/// 
/// # Examples
/// ```
/// use sdomain_test::sdomain::{Fs, parallel};
/// 
/// let model1 = Fs{
///     numer:vec![1.0, 2.0, 3.0],
///     denom:vec![-4.0, 5.0, 2.0, -1.0],
///     prec:1
/// };
/// let model2 = Fs{
///     numer:vec![3.0, -5.0],
///     denom:vec![-3.0],
///     prec:1
/// };
/// 
/// let result = parallel(model1, model2);
/// 
/// assert_eq!("( 3.0e0 + 1.0e0s + -1.0e0s² + -1.5e1s³ ) / ( -1.5e1 + 2.9e1s + -2.8e1s² + -1.3e1s³ + 5.0e0s⁴ )", format!("{}", result));
/// ```
pub fn parallel(fs1: Fs, fs2: Fs) -> Fs {
    (fs1.invert() + &fs2.invert()).invert()
}

/// Combine a collection of models assuming they are all in parallel.
/// Returns a combined model.
/// Returns an `Option<Fs>`, where None means the collection was empty.
/// ```text
/// 
///           ┌───┴───┬──────┬──...───┐  
///        ┌──┴─┐  ┌──┴─┐ ┌──┴─┐   ┌──┴─┐
/// T(s) = │F(s)│  │G(s)│ │H(s)│...│Z(s)│
///        └──┬─┘  └──┬─┘ └──┬─┘   └──┬─┘
///           └───┬───┴──────┴──...───┘  
/// 
/// ```
/// 
/// # Arguments
/// * `iter` - An iterator pointing to the models to combine.
/// 
/// # Examples
/// ```
/// use sdomain_test::sdomain::{Fs, parallel_multi};
/// 
/// let models = [
///     Fs{
///         numer:vec![1.0, 2.0, 3.0],
///         denom:vec![-4.0, 5.0, 2.0, -1.0],
///         prec:1
///     },
///     Fs{
///         numer:vec![0.0, 4.5],
///         denom:vec![1.0],
///         prec:1
///     },
///     Fs{
///         numer:vec![1.0],
///         denom:vec![2.5, -7.0],
///         prec:1
///     },
///     Fs{
///         numer:vec![1.0],
///         denom:vec![1.0],
///         prec:1
///     }
/// ];
/// 
/// let result = parallel_multi(models.into_iter());
/// 
/// assert_eq!("( 0.0e0 + 4.5e0s + 9.0e0s² + 1.4e1s³ ) / ( 1.0e0 + -2.5e-1s + 2.6e1s² + -6.8e0s³ + -9.9e1s⁴ )", format!("{}", result));
/// ```
pub fn parallel_multi<T: Iterator<Item = Fs> + Clone>(iter: T) -> Fs {
    // if iter.clone().peekable().peek().is_none() {return None;}
    // Some(iter.map(|fs| fs.invert()).sum::<Fs>().invert())
    iter.map(|fs| fs.invert()).sum::<Fs>().invert()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn init_test() {
        let fs = Fs{numer:vec![1.0, 0.0, 2.0, 3.0], denom:vec![5.0, 6.0, 7.0, 8.0], prec: 4};
        assert_eq!(
            "[1.0, 0.0, 2.0, 3.0] [5.0, 6.0, 7.0, 8.0] 4",
            format!("{0:?} {1:?} {2}", fs.numer, fs.denom, fs.prec)
        );
    }

    #[test]
    fn merge_rad_test() {
        //                                               1     2    4    8    16    32    64  128  256
        let c = Fs::merge_rad(&vec![1.0, 2.0, -3.0, 4.0, 5.0, -6.0, -7.0, 8.0, 9.0], 2.0, 4);
        assert_eq!("2.8450e3 + -1.2440e3i", format!("{c}"));
    }

    #[test]
    fn merge_freq_test() {
        //                                                1     2    4    8    16    32    64  128  256
        let c = Fs::merge_freq(&vec![1.0, 2.0, -3.0, 4.0, 5.0, -6.0, -7.0, 8.0, 9.0], 1.0/PI, 4);
        assert_eq!("2.8450e3 + -1.2440e3i", format!("{c}"));
    }

    #[test]
    fn calculate_rad_test() {
        let converted = Fs{
                numer:vec![1.0, -2.0, 3.0, 4.0, 5.0],                  //  +3.0 + -6.0i
                denom:vec![-1.0, 6.0, 0.5, 4.0, -8.0, -2.0, 4.0, 5.0], // -13.5 + -5.0i (207.25)
                prec :3
        }.calculate_rad(1.0); // ( -10.5 + 96i ) / 207.25
        assert_eq!("-5.066e-2 + 4.632e-1i", format!("{converted}"));
    }

    #[test]
    fn calculate_freq_test() {
        let converted = Fs{
                numer:vec![1.0, -2.0, 3.0, 4.0, 5.0],                  //  +3.0 + -6.0i
                denom:vec![-1.0, 6.0, 0.5, 4.0, -8.0, -2.0, 4.0, 5.0], // -13.5 + -5.0i (207.25)
                prec :3
        }.calculate_freq(0.5/PI); // ( -10.5 + 96i ) / 207.25
        assert_eq!("-5.066e-2 + 4.632e-1i", format!("{converted}"));
    }

    #[test]
    fn scale_test() {
        let scaled = Fs{
            numer:vec![1.0, -2.0],
            denom:vec![-1.0],
            prec :1
        }.scale(4.0);
        assert_eq!("( 4.0e0 + -8.0e0s ) / ( -1.0e0 )", format!("{scaled}"));
    }

    #[test]
    fn invert_test() {
        let inverted = Fs{
            numer:vec![1.0, -2.0],
            denom:vec![-1.0],
            prec :1
        }.invert();
        assert_eq!("( -1.0e0 ) / ( 1.0e0 + -2.0e0s )", format!("{inverted}"));
    }

    #[test]
    fn trim_helper_test() {
        let mut polynomial = vec![1.0, 0.0, 2.0, 0.0, 0.0];
        Fs::trim_helper(&mut polynomial);
        assert_eq!("[1.0, 0.0, 2.0]", format!("{:?}", polynomial));
        let mut collection2 = vec![0.0, 0.0];
        Fs::trim_helper(&mut collection2);
        assert_eq!("[0.0]", format!("{:?}", collection2));
    }

    #[test]
    fn trim_leading_zeroes_test() {
        let mut fs = Fs{
            numer: vec![0.0, 0.0, 0.0, 0.0],
            denom: vec![2.0, 3.0, 0.0, 0.0, 4.0, 0.0],
            prec: 1
        };
        fs.trim_leading_zeroes();
        assert_eq!("( 0.0e0 ) / ( 2.0e0 + 3.0e0s + 0.0e0s² + 0.0e0s³ + 4.0e0s⁴ )", format!("{}", fs));
    }
}

#[cfg(test)]
mod conversion_tests {
    use super::*;
    
    #[test]
    fn fs2_complex_rad_test() {
        let converted = Complex::from_rad(
            Fs{
                numer:vec![1.0, -2.0, 3.0, 4.0, 5.0],                  //  +3.0 + -6.0i
                denom:vec![-1.0, 6.0, 0.5, 4.0, -8.0, -2.0, 4.0, 5.0], // -13.5 + -5.0i (207.25)
                prec :3
            },
            1.0
        ); // ( -10.5 + 96i ) / 207.25
        assert_eq!("-5.066e-2 + 4.632e-1i", format!("{converted}"));
    }

    #[test]
    fn fs2_complex_freq_test() {
        let converted = Complex::from_freq(
            Fs{
                numer:vec![1.0, -2.0, 3.0, 4.0, 5.0],                  //  +3.0 + -6.0i
                denom:vec![-1.0, 6.0, 0.5, 4.0, -8.0, -2.0, 4.0, 5.0], // -13.5 + -5.0i (207.25)
                prec :3
            },
            0.5/PI
        ); // ( -10.5 + 96i ) / 207.25
        assert_eq!("-5.066e-2 + 4.632e-1i", format!("{converted}"));
    }
}


#[cfg(test)]
mod display_tests {
    use super::*;

    #[test]
    fn fmt_s_term_0_test() {
        let res = Fs::fmt_s_term(0);
        assert_eq!("", res);
    }

    #[test]
    fn fmt_s_term_1_test() {
        let res = Fs::fmt_s_term(1);
        assert_eq!("s", res);
    }

    #[test]
    fn fmt_s_term_var_test() {
        let res = Fs::fmt_s_term(45);
        assert_eq!("s⁴⁵", res);
    }

    #[test]
    fn to_str_test() {
        let fs = Fs{
            numer:vec![1.0, -2.0, 3.0, 4.0, 5.0],
            denom:vec![-1.0, 6.0, 0.5, 4.0, -8.0, -2.0, 4.0, 5.0],
            prec :0
        };
        assert_eq!("1e0 + -2e0s + 3e0s² + 4e0s³ + 5e0s⁴", fs.to_str_numer(false));
        assert_eq!("1e0s⁰ + -2e0s¹ + 3e0s² + 4e0s³ + 5e0s⁴", fs.to_str_numer(true));
        assert_eq!("-1e0 + 6e0s + 5e-1s² + 4e0s³ + -8e0s⁴ + -2e0s⁵ + 4e0s⁶ + 5e0s⁷", fs.to_str_denom(false));
        assert_eq!("-1e0s⁰ + 6e0s¹ + 5e-1s² + 4e0s³ + -8e0s⁴ + -2e0s⁵ + 4e0s⁶ + 5e0s⁷", fs.to_str_denom(true));
    }

    #[test]
    fn fmt_test() {
        let fs = Fs{
            numer:vec![1.0, -2.0, 3.0, 4.0, 5.0],
            denom:vec![-1.0, 6.0, 0.5, 4.0, -8.0, -2.0, 4.0, 5.0],
            prec :0
        };
        assert_eq!("( 1e0 + -2e0s + 3e0s² + 4e0s³ + 5e0s⁴ ) / ( -1e0 + 6e0s + 5e-1s² + 4e0s³ + -8e0s⁴ + -2e0s⁵ + 4e0s⁶ + 5e0s⁷ )", format!("{}", fs));
    }
}


#[cfg(test)]
mod arith_tests {
    use super::*;

    fn gen_factors() -> (Fs, Fs) {
        (
            Fs{numer:vec![1.0, 2.0], denom:vec![4.0], prec:4},
            Fs{numer:vec![2.5, -0.25, 0.25], denom:vec![3.0, 0.5], prec:1}
        )
    }

    #[test]
    fn neg_test() {
        let (lhs, _) = gen_factors();
        assert_eq!(
            "( -1.0000e0 + -2.0000e0s ) / ( 4.0000e0 )",
            format!("{}", -lhs)
        );
    }

    #[test]
    fn accumulate_test() {
        let bucket = vec![
            vec![],
            vec![0.5, -3.0],
            vec![1.0, 2.0, 3.0, 4.0],
            vec![-1.5]
        ];
        assert_eq!("[0.0, -1.0, 3.0, 4.0]", format!("{:?}", Fs::accumulate(bucket)));
    }

    #[test]
    fn multiply_helper_test() {
        let lhs = vec![0.5, 3.0];
        let rhs = vec![4.0, -2.0, 1.0];
        // 2.0  -1.0    0.5
        // 0.0  12.0    -6.0    3.0
        // ------------------------
        // 2.0  11.0    -5.5    3.0
        assert_eq!("[2.0, 11.0, -5.5, 3.0]", format!("{:?}", Fs::multiply_helper(&lhs, &rhs)));
    }

    #[test]
    fn mul_test() {
        let (lhs, rhs) = gen_factors();
        // numer
        // 2.50 -0.25   0.25
        // 0.00 5.00    -0.50   0.5
        // ------------------------
        // 2.50 4.75    -0.25   0.5
        //denom
        // 12.0 2.0
        assert_eq!(
            "( 2.5e0 + 4.8e0s + -2.5e-1s² + 5.0e-1s³ ) / ( 1.2e1 + 2.0e0s )",
            format!("{}", lhs * &rhs)
        );
    }

    #[test]
    fn div_test() {
        let (lhs, rhs) = gen_factors();
        // numer
        // 3.0  0.5
        // 0.0  6.0 1.0
        // ------------
        // 3.0  6.5 1.0
        //denom
        // 10.0 -1.0 1.0
        assert_eq!(
            "( 3.0e0 + 6.5e0s + 1.0e0s² ) / ( 1.0e1 + -1.0e0s + 1.0e0s² )",
            format!("{}", lhs / &rhs)
        );
    }

    #[test]
    fn add_eq_denom_test() {
        let (lhs, mut rhs) = gen_factors();
        rhs.denom = lhs.denom.clone();

        assert_eq!(
            "( 3.5e0 + 1.8e0s + 2.5e-1s² ) / ( 4.0e0 )",
            format!("{}", lhs + &rhs)
        );
    }

    #[test]
    fn add_diff_denom_test() {
        let (lhs, rhs) = gen_factors();

        // lhs.numer*rhs.denom
        // 3.0  0.5
        // 0.0  6.0 1.0
        // ------------
        // 3.0  6.5 1.0
        // lhs.denom*rhs.numer
        // 10.0 -1.0 1.0
        // lhs.numer*rhs.denom + lhs.denom*rhs.numer
        // 13.0 5.5 2.0
        // lhs.denom*rhs.denom
        // 12.0 2.0
        assert_eq!(
            "( 1.3e1 + 5.5e0s + 2.0e0s² ) / ( 1.2e1 + 2.0e0s )",
            format!("{}", lhs + &rhs)
        );
    }

    #[test]
    fn sub_test() {
        let (lhs, rhs) = gen_factors();

        // lhs.numer*rhs.denom + lhs.denom*rhs.numer
        // 3.0      6.5 1.0
        // -10.0    1.0 -1.0
        // -----------------
        // -7.0     7.5 0.0
        // lhs.denom*rhs.denom
        // 12.0 2.0
        assert_eq!(
            "( -7.0e0 + 7.5e0s ) / ( 1.2e1 + 2.0e0s )",
            format!("{}", lhs - &rhs)
        );
    }

    #[test]
    fn mul_assign_test() {
        let (mut lhs, rhs) = gen_factors();
        lhs *= &rhs;
        // numer
        // 2.50 -0.25   0.25
        // 0.00 5.00    -0.50   0.5
        // ------------------------
        // 2.50 4.75    -0.25   0.5
        //denom
        // 12.0 2.0
        assert_eq!(
            "( 2.5e0 + 4.8e0s + -2.5e-1s² + 5.0e-1s³ ) / ( 1.2e1 + 2.0e0s )",
            format!("{}", lhs)
        );
    }

    #[test]
    fn div_assign_test() {
        let (mut lhs, rhs) = gen_factors();
        lhs /= &rhs;
        // numer
        // 3.0  0.5
        // 0.0  6.0 1.0
        // ------------
        // 3.0  6.5 1.0
        //denom
        // 10.0 -1.0 1.0
        assert_eq!(
            "( 3.0e0 + 6.5e0s + 1.0e0s² ) / ( 1.0e1 + -1.0e0s + 1.0e0s² )",
            format!("{}", lhs)
        );
    }

    #[test]
    fn add_assign_eq_denom_test() {
        let (mut lhs, mut rhs) = gen_factors();
        rhs.denom = lhs.denom.clone();
        lhs += &rhs;

        assert_eq!(
            "( 3.5e0 + 1.8e0s + 2.5e-1s² ) / ( 4.0e0 )",
            format!("{}", lhs)
        );
    }

    #[test]
    fn add_assign_diff_denom_test() {
        let (mut lhs, rhs) = gen_factors();
        lhs += &rhs;

        // lhs.numer*rhs.denom
        // 3.0  0.5
        // 0.0  6.0 1.0
        // ------------
        // 3.0  6.5 1.0
        // lhs.denom*rhs.numer
        // 10.0 -1.0 1.0
        // lhs.numer*rhs.denom + lhs.denom*rhs.numer
        // 13.0 5.5 2.0
        // lhs.denom*rhs.denom
        // 12.0 2.0
        assert_eq!(
            "( 1.3e1 + 5.5e0s + 2.0e0s² ) / ( 1.2e1 + 2.0e0s )",
            format!("{}", lhs)
        );
    }

    #[test]
    fn sub_assign_test() {
        let (mut lhs, rhs) = gen_factors();
        lhs -= &rhs;

        // lhs.numer*rhs.denom + lhs.denom*rhs.numer
        // 3.0      6.5 1.0
        // -10.0    1.0 -1.0
        // -----------------
        // -7.0     7.5 0.0
        // lhs.denom*rhs.denom
        // 12.0 2.0
        assert_eq!(
            "( -7.0e0 + 7.5e0s ) / ( 1.2e1 + 2.0e0s )",
            format!("{}", lhs)
        );
    }
}


#[cfg(test)]
mod gen_tests {
    use super::*;

    #[test]
    fn gen_s_test() {
        assert_eq!("( 0.0e0 + 1.0e0s ) / ( 1.0e0 )", format!("{}", gen::s_prec(1)));
    }

    #[test]
    fn gen_step_test() {
        assert_eq!("( 1.0e0 ) / ( 0.0e0 + 1.0e0s )", format!("{}", gen::step_prec(1)));
    }

    #[test]
    fn gen_resistor_test() {
        assert_eq!("( 4.36e1 ) / ( 1.00e0 )", format!("{}", gen::resistor_prec(43.6, 2)));
    }

    #[test]
    fn gen_capacitor_test() {
        assert_eq!("( 1.00e0 ) / ( 0.00e0 + 4.36e1s )", format!("{}", gen::capacitor_prec(43.6, 2)));
    }

    #[test]
    fn gen_inductor_test() {
        assert_eq!("( 0.00e0 + 4.36e1s ) / ( 1.00e0 )", format!("{}", gen::inductor_prec(43.6, 2)));
    }

    #[test]
    fn gen_rcl_test() {
        assert_eq!("( 1.0e0 + 2.2e1s + 6.0e0s² ) / ( 0.0e0 + 2.0e0s )", format!("{}", gen::rcl_prec(11.0, 2.0, 3.0, 1)));
    }
}

#[cfg(test)]
mod parallel_tests {
    use super::*;

    #[test]
    fn parallel_test() {
        let model1 = Fs{
            numer:vec![1.0],
            denom:vec![-4.0, 1.0],
            prec:1
        };
        let model2 = Fs{
            numer:vec![3.0, -5.0],
            denom:vec![-3.0],
            prec:1
        };
        // denom
        // -3 - 12 + 20s + 3s - 5s²
        // -15 + 23s - 5s²
        // numer
        // 3 - 5s
        let result = parallel(model1, model2);
        
        assert_eq!("( 3.0e0 + -5.0e0s ) / ( -1.5e1 + 2.3e1s + -5.0e0s² )", format!("{}", result));
    }

    #[test]
    fn parallel_short_test() {
        let lhs = Fs{
            numer:vec![1.0],
            denom:vec![-4.0, 1.0],
            prec:1
        };
        let result = parallel(lhs.clone(), gen::zero_prec(100));
        assert_eq!("( 0.0e0 ) / ( 1.0e0 )", format!("{}", result));
        let result2 = parallel(gen::zero_prec(100), lhs);
        assert_eq!("( 0.0e0 ) / ( 1.0e0 )", format!("{}", result2));
    }

    #[test]
    fn parallel_multishort_test() {
        let result = parallel(gen::zero_prec(1), gen::zero_prec(1));
        assert_eq!("( 0.0e0 ) / ( 2.0e0 )", format!("{}", result));
        let result2 = parallel(Fs{numer:vec![0.0], denom:vec![2.0, 1.0], prec:1}, gen::zero_prec(1));
        assert_eq!("( 0.0e0 ) / ( 3.0e0 + 1.0e0s )", format!("{}", result2));
        let mut lhs = Fs{numer:vec![0.0, 0.0], denom:vec![2.0, 1.0], prec:1};
        lhs.trim_leading_zeroes();
        let result2 = parallel(lhs, gen::zero_prec(1));
        assert_eq!("( 0.0e0 ) / ( 3.0e0 + 1.0e0s )", format!("{}", result2));
    }

    #[test]
    fn parallel_opencircuit_test() {
        let lhs = Fs{
            numer:vec![1.0],
            denom:vec![-4.0, 1.0],
            prec:1
        };
        let result = parallel(lhs.clone(), gen::zero_prec(100).invert());
        
        assert_eq!(format!("{}", lhs), format!("{}", result));
    }

    #[test]
    fn parallel_multi_test() {
        let models = [
            Fs{
                numer:vec![1.0],
                denom:vec![-4.0, 1.0],
                prec:1
            },
            Fs{
                numer:vec![3.0, -5.0],
                denom:vec![-3.0],
                prec:1
            },
            Fs{
                numer:vec![4.0],
                denom:vec![1.0],
                prec:1
            },
        ];
        // First sum
        //  numer
        //  -15 + 23s - 5s²
        //  denom
        //  3 - 5s
        // Second sum
        //  numer
        //  -60 + 92s - 20s² + 3 - 5s
        //  -57 + 87s - 20s²
        //  denom
        //  12 - 20s
        // Invert
        let result = parallel_multi(models.into_iter());
        
        assert_eq!("( 1.2e1 + -2.0e1s ) / ( -5.7e1 + 8.7e1s + -2.0e1s² )", format!("{}", result));
    }

    #[test]
    fn parallel_multi_shortcircuit_test() {
        let models = [
            Fs{
                numer:vec![1.0],
                denom:vec![-4.0, 1.0],
                prec:1
            },
            gen::zero_prec(100),
            Fs{
                numer:vec![3.0, -5.0],
                denom:vec![-3.0],
                prec:1
            },
            gen::zero_prec(100),
            Fs{
                numer:vec![4.0],
                denom:vec![1.0],
                prec:1
            },
        ];
        let result = parallel_multi(models.into_iter());
        //( 3.0e0 + -5.0e0s + 0.0e0s² ) / ( 0.0e0 + 0.0e0s )
        //( 1.0e0 ) / ( 0.0e0 )
        assert_eq!("( 0.0e0 ) / ( 1.6e1 + -2.0e1s )", format!("{}", result));
    }

    #[test]
    fn parallel_multi_opencircuit_test() {
        let models = [
            Fs{
                numer:vec![1.0],
                denom:vec![-4.0, 1.0],
                prec:1
            },
            gen::zero_prec(100).invert(),
            Fs{
                numer:vec![3.0, -5.0],
                denom:vec![-3.0],
                prec:1
            },
            gen::zero_prec(100).invert(),
            Fs{
                numer:vec![4.0],
                denom:vec![1.0],
                prec:1
            },
        ];
        let result = parallel_multi(models.into_iter());
        
        assert_eq!("( 1.2e1 + -2.0e1s ) / ( -5.7e1 + 8.7e1s + -2.0e1s² )", format!("{}", result));
    }
}