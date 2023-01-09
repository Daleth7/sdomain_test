use std::{fmt::Display, f64::consts::PI, ops::{Neg, Add, Sub, Mul, Div}, cmp::min};

use crate::complex::Complex;

#[derive(Debug, PartialEq)]
pub struct Fs {
    pub numer: Vec<f64>,
    pub denom: Vec<f64>,
    pub prec: usize,
}

impl Fs {
    pub fn merge_rad(collection: &Vec<f64>, rad: f64, precision: usize) -> Complex {
        let mut real = 0.0;
        let mut imag = 0.0;
        let mut val = 1.0;
        for (idx, coeff) in collection.iter().enumerate() {
            if idx % 2 == 0 { real += coeff*val*(if idx % 4 == 0 {1.0} else {-1.0}); }
            else            { imag += coeff*val*(if (idx-1) % 4 == 0 {1.0} else {-1.0}); }
            val *= rad;
        }
        Complex { real, imag, prec: precision }
    }

    pub fn merge_freq(collection: &Vec<f64>, freq: f64, precision: usize) -> Complex {
        Self::merge_rad(collection, freq*2.0*PI, precision)
    }

    fn fmt_s_term(exponent: usize) -> String {
        match exponent {
            0 => String::from(""),
            1 => String::from("s"),
            _ => format!("s^{exponent}"),
        }
    }

    pub fn to_str(collection: &Vec<f64>, precision: usize, show_all: bool) -> String {
        let terms: Vec<String> = collection.iter().enumerate().map(|(idx, val)| {
            if show_all {
                format!("{val:.0$e}s^{1}", precision, idx)
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

    pub fn to_str_numer(&self, show_all: bool) -> String {
        Self::to_str(&self.numer, self.prec, show_all)
    }

    pub fn to_str_denom(&self, show_all: bool) -> String {
        Self::to_str(&self.denom, self.prec, show_all)
    }

    pub fn calculate_rad(&self, rad: f64) -> Complex {
        Fs::merge_rad(&self.numer, rad, self.prec) / Fs::merge_rad(&self.denom, rad, self.prec)
    }

    pub fn calculate_freq(&self, freq: f64) -> Complex {
        Fs::merge_freq(&self.numer, freq, self.prec) / Fs::merge_freq(&self.denom, freq, self.prec)
    }

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

    fn multiply_helper(lhs: Vec<f64>, rhs: Vec<f64>) -> Vec<f64> {
        let mut bucket: Vec::<Vec::<f64>> = rhs.into_iter()
            .map(|rhs_coeff| {
                lhs.clone()
                    .into_iter()
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
}

impl Display for Fs {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let numer_str = self.to_str_numer(false);
        let denom_str = self.to_str_denom(false);
        write!(f, "( {0} ) / ( {1} )", numer_str, denom_str)
    }
}

impl Complex {
    fn from_rad(value: Fs, rad: f64) -> Self {
        value.calculate_rad(rad)
    }

    fn from_freq(value: Fs, freq: f64) -> Self {
        value.calculate_freq(freq)
    }
}

impl Neg for Fs {
    type Output = Fs;

    fn neg(self) -> Self::Output {
        Self{numer:self.numer.iter().map(|x| x*-1.0).collect(), denom:self.denom, prec:self.prec}
    }
}

impl Add for Fs {
    type Output = Fs;

    fn add(self, rhs: Self) -> Self::Output {
        if self.denom == rhs.denom {
            Self{
                numer: Self::accumulate([self.numer, rhs.numer]),
                denom: self.denom,
                prec: min(self.prec, rhs.prec)
            }
        } else {
            let lhs_denom = self.denom.clone();
            let rhs_denom = rhs.denom.clone();
            Self{
                denom: Self::multiply_helper(self.denom , rhs.denom),
                numer: Self::accumulate([
                    Self::multiply_helper(self.numer, rhs_denom),
                    Self::multiply_helper(rhs.numer, lhs_denom)
                    ]),
                prec: min(self.prec, rhs.prec)
            }
        }
    }
}

impl Sub for Fs {
    type Output = Fs;

    fn sub(self, rhs: Self) -> Self::Output {
        self + (-rhs)
    }
}

impl Mul for Fs {
    type Output = Fs;

    fn mul(self, rhs: Self) -> Self::Output {
        Self{
            numer: Self::multiply_helper(self.numer, rhs.numer),
            denom: Self::multiply_helper(self.denom, rhs.denom),
            prec: min(self.prec, rhs.prec),
        }
    }
}

impl Div for Fs {
    type Output = Fs;

    fn div(self, rhs: Self) -> Self::Output {
        Self{
            numer: Self::multiply_helper(self.numer, rhs.denom),
            denom: Self::multiply_helper(self.denom, rhs.numer),
            prec: min(self.prec, rhs.prec),
        }
    }
}

pub mod gen {
    use super::Fs;

    pub fn s(prec: usize) -> Fs {Fs{numer:vec![0.0, 1.0], denom:vec![1.0], prec}}
    pub fn step(prec: usize) -> Fs {Fs{numer:vec![1.0], denom:vec![0.0, 1.0], prec}}
    pub fn resistor(r: f64, prec: usize) -> Fs {Fs{numer:vec![r], denom:vec![1.0], prec}}
    pub fn capacitor(c: f64, prec: usize) -> Fs {Fs{numer:vec![1.0], denom:vec![0.0, c], prec}}
    pub fn inductor(l: f64, prec: usize) -> Fs {Fs{numer:vec![0.0, l], denom:vec![1.0], prec}}

    pub fn rcl(r: f64, c: f64, l: f64, prec: usize) -> Fs {Fs{numer:vec![1.0, r*c, l*c], denom:vec![0.0, c], prec}}
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
        assert_eq!("2845.0000 + -1244.0000i", format!("{c}"));
    }

    #[test]
    fn merge_freq_test() {
        //                                                1     2    4    8    16    32    64  128  256
        let c = Fs::merge_freq(&vec![1.0, 2.0, -3.0, 4.0, 5.0, -6.0, -7.0, 8.0, 9.0], 1.0/PI, 4);
        assert_eq!("2845.0000 + -1244.0000i", format!("{c}"));
    }

    #[test]
    fn calculate_rad_test() {
        let converted = Fs{
                numer:vec![1.0, -2.0, 3.0, 4.0, 5.0],                  //  +3.0 + -6.0i
                denom:vec![-1.0, 6.0, 0.5, 4.0, -8.0, -2.0, 4.0, 5.0], // -13.5 + -5.0i (207.25)
                prec :3
        }.calculate_rad(1.0); // ( -10.5 + 96i ) / 207.25
        assert_eq!("-0.051 + 0.463i", format!("{converted}"));
    }

    #[test]
    fn calculate_freq_test() {
        let converted = Fs{
                numer:vec![1.0, -2.0, 3.0, 4.0, 5.0],                  //  +3.0 + -6.0i
                denom:vec![-1.0, 6.0, 0.5, 4.0, -8.0, -2.0, 4.0, 5.0], // -13.5 + -5.0i (207.25)
                prec :3
        }.calculate_freq(0.5/PI); // ( -10.5 + 96i ) / 207.25
        assert_eq!("-0.051 + 0.463i", format!("{converted}"));
    }

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
        assert_eq!("-0.051 + 0.463i", format!("{converted}"));
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
        assert_eq!("-0.051 + 0.463i", format!("{converted}"));
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
        assert_eq!("s^45", res);
    }

    #[test]
    fn to_str_test() {
        let fs = Fs{
            numer:vec![1.0, -2.0, 3.0, 4.0, 5.0],
            denom:vec![-1.0, 6.0, 0.5, 4.0, -8.0, -2.0, 4.0, 5.0],
            prec :0
        };
        assert_eq!("1e0 + -2e0s + 3e0s^2 + 4e0s^3 + 5e0s^4", fs.to_str_numer(false));
        assert_eq!("1e0s^0 + -2e0s^1 + 3e0s^2 + 4e0s^3 + 5e0s^4", fs.to_str_numer(true));
        assert_eq!("-1e0 + 6e0s + 5e-1s^2 + 4e0s^3 + -8e0s^4 + -2e0s^5 + 4e0s^6 + 5e0s^7", fs.to_str_denom(false));
        assert_eq!("-1e0s^0 + 6e0s^1 + 5e-1s^2 + 4e0s^3 + -8e0s^4 + -2e0s^5 + 4e0s^6 + 5e0s^7", fs.to_str_denom(true));
    }

    #[test]
    fn fmt_test() {
        let fs = Fs{
            numer:vec![1.0, -2.0, 3.0, 4.0, 5.0],
            denom:vec![-1.0, 6.0, 0.5, 4.0, -8.0, -2.0, 4.0, 5.0],
            prec :0
        };
        assert_eq!("( 1e0 + -2e0s + 3e0s^2 + 4e0s^3 + 5e0s^4 ) / ( -1e0 + 6e0s + 5e-1s^2 + 4e0s^3 + -8e0s^4 + -2e0s^5 + 4e0s^6 + 5e0s^7 )", format!("{}", fs));
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
        assert_eq!("[2.0, 11.0, -5.5, 3.0]", format!("{:?}", Fs::multiply_helper(lhs, rhs)));
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
            "( 2.5e0 + 4.8e0s + -2.5e-1s^2 + 5.0e-1s^3 ) / ( 1.2e1 + 2.0e0s )",
            format!("{}", lhs * rhs)
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
            "( 3.0e0 + 6.5e0s + 1.0e0s^2 ) / ( 1.0e1 + -1.0e0s + 1.0e0s^2 )",
            format!("{}", lhs / rhs)
        );
    }

    #[test]
    fn add_eq_denom_test() {
        let (lhs, mut rhs) = gen_factors();
        rhs.denom = lhs.denom.clone();

        assert_eq!(
            "( 3.5e0 + 1.8e0s + 2.5e-1s^2 ) / ( 4.0e0 )",
            format!("{}", lhs + rhs)
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
            "( 1.3e1 + 5.5e0s + 2.0e0s^2 ) / ( 1.2e1 + 2.0e0s )",
            format!("{}", lhs + rhs)
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
            "( -7.0e0 + 7.5e0s + 0.0e0s^2 ) / ( 1.2e1 + 2.0e0s )",
            format!("{}", lhs - rhs)
        );
    }
}


#[cfg(test)]
mod gen_tests {
    use super::*;

    #[test]
    fn gen_s_test() {
        assert_eq!("( 0.0e0 + 1.0e0s ) / ( 1.0e0 )", format!("{}", gen::s(1)));
    }

    #[test]
    fn gen_step_test() {
        assert_eq!("( 1.0e0 ) / ( 0.0e0 + 1.0e0s )", format!("{}", gen::step(1)));
    }

    #[test]
    fn gen_resistor_test() {
        assert_eq!("( 4.36e1 ) / ( 1.00e0 )", format!("{}", gen::resistor(43.6, 2)));
    }

    #[test]
    fn gen_capacitor_test() {
        assert_eq!("( 1.00e0 ) / ( 0.00e0 + 4.36e1s )", format!("{}", gen::capacitor(43.6, 2)));
    }

    #[test]
    fn gen_inductor_test() {
        assert_eq!("( 0.00e0 + 4.36e1s ) / ( 1.00e0 )", format!("{}", gen::inductor(43.6, 2)));
    }

    #[test]
    fn gen_rcl_test() {
        assert_eq!("( 1.0e0 + 2.2e1s + 6.0e0s^2 ) / ( 0.0e0 + 2.0e0s )", format!("{}", gen::rcl(11.0, 2.0, 3.0, 1)));
    }
}