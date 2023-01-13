use std::f64::consts::PI;

use crate::passives::packages::Package;
use crate::sdomain::{self, Fs};

#[derive(Debug, Copy, Clone, PartialEq)]
pub struct Capacitor {
    pub val: f64,
    pub esr: f64,
    pub esl: f64,
    pub lmnt: f64,
    pub package: Package,
}

impl Capacitor {
    pub fn resonant(&self) -> f64 {
        1.0/(self.val*(self.esl+self.lmnt)).sqrt()/2.0/PI
    }

    pub fn model(&self) -> Fs {
        sdomain::gen::rcl(self.esr, self.val, self.esl+self.lmnt)
    }

    pub fn from(val: f64, package: &str) -> Capacitor {
        let database = include!("capacitor_database.in");
        let package_convert = match package.parse() {
            Ok(val) => val,
            Err(_) => panic!("Unknown package specified ({package})!")
        };
        let iter = database
            .iter()
            .filter(|c|
                c.package == package_convert
             && val < (c.val+1e-12)
             && val > (c.val-1e-12)
            ).next();
        match iter {
            Some(cap) => *cap,
            None => panic!("Could not find capacitor with {val:.3e}F capacitance in {package} package!")
        }
    }

    pub fn from_resonant(center: f64, err: f64) -> Option<Capacitor> {
        let database = include!("capacitor_database.in");
        let result = database
            .iter()
            .filter(|c| {
                let resonant = c.resonant();
                resonant < (center+err) && resonant > (center-err)
            }).next();
        Some(*result?)
    }
}

impl std::fmt::Display for Capacitor {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "Capacitor {{ val: {0:.3}nF, esr: {1:.1}mΩ, esl: {2:.1}nH, lmnt: {3:.1}nH, package: {4:?} }}",
            self.val*1e9,
            self.esr*1e3,
            self.esl*1e9,
            self.lmnt*1e9,
            self.package
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn resonant_test() {
        let cap = Capacitor{val:3.0e-6, esr:7.0, esl:2.0e-9, lmnt:17.0e-9, package:Package::Electrolytic};
        let resonant = cap.resonant();
        let expected = 666626.0;
        let err = 1.0;
        println!("Resonant = {resonant:.1}Hz | expected = {expected:.1}Hz");
        assert!(resonant < (expected+err) && resonant > (expected-err));
    }

    #[test]
    fn model_test() {
        let cap = Capacitor{val:1.0, esr:2.0, esl:3.0, lmnt:4.0, package:Package::Electrolytic};
        assert_eq!("( 1.000e0 + 2.000e0s + 7.000e0s^2 ) / ( 0.000e0 + 1.000e0s )", format!("{}", cap.model()))
    }

    #[test]
    fn from_test() {
        let cap1 = Capacitor::from(4.7e-6, "0805");
        assert_eq!("( 1.000e0 + 1.880e-8s + 7.050e-15s^2 ) / ( 0.000e0 + 4.700e-6s )", format!("{}", cap1.model()));
        let cap2 = Capacitor::from(4.7e-6, "2012M");
        assert_eq!("( 1.000e0 + 1.880e-8s + 7.050e-15s^2 ) / ( 0.000e0 + 4.700e-6s )", format!("{}", cap2.model()));
    }

    #[test]
    fn from_resonant_test() {
        let cap1 = Capacitor::from_resonant(28.677e6, 0.01e6).unwrap();
        assert_eq!("( 1.000e0 + 9.460e-10s + 3.080e-17s^2 ) / ( 0.000e0 + 2.200e-8s )", format!("{}", cap1.model()));
        let cap2 = Capacitor::from_resonant(28.677e6, 1.0);
        assert_eq!(None, cap2);
    }

    #[test]
    #[should_panic]
    fn from_panic_package_test() {
        Capacitor::from(0.5e-6, "Unknown");
    }

    #[test]
    #[should_panic]
    fn from_panic_val_test() {
        Capacitor::from(0.5e-6, "0805");
    }

    #[test]
    fn display_test() {
        let cap = Capacitor{val:3.0e-6, esr:0.123, esl:2.0e-9, lmnt:1.0e-9, package:Package::SMD_4532M};
        assert_eq!("Capacitor { val: 3000.000nF, esr: 123.0mΩ, esl: 2.0nH, lmnt: 1.0nH, package: SMD_4532M }", format!("{}", cap));
    }
}