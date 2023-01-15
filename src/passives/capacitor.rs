use std::f64::consts::PI;

use crate::passives::packages::Package;
use crate::sdomain::{self, Fs};

/// Represents the capacitor electrical component. Includes paraistic parameters such as ESR and ESL.
/// 
/// # Examples
/// ```
/// use complextest::passives::capacitor::Capacitor;
/// use complextest::passives::packages::Package;
/// 
/// let mut cap = Capacitor{
///     val: 10e-6,   //F
///     esr: 0.040,   //Ω
///     esl: 0.4e-9,  //H
///     lmnt: 1.0e-9, //H
///     package: Package::SMD_2012M
/// };
/// 
/// // Prints "Capacitor { val: 10000.000nF, esr: 40.0mΩ, esl: 0.4nH, lmnt: 1.0nH, package: SMD_2012M }"
/// println!("{cap}");
/// 
/// // Use a template with a certain capacitance and case code. The ESR and ESL default come from the template.
/// let cap2 = Capacitor::from(22e-9, "0201"); // Or "0603M" for metric
/// // Prints "Capacitor { val: 22.000nF, esr: 39.0mΩ, esl: 0.3nH, lmnt: 1.0nH, package: SMD_0603M }"
/// println!("{cap2}");
/// 
/// // Generate an impedance model in the s-domain from the cap. parameters
/// // Prints "( 1.000e0 + 8.580e-10s + 2.860e-17s² ) / ( 0.000e0 + 2.200e-8s )"
/// println!("{}", cap2.model());
/// 
/// // Update capacitance value with derated value from manufacturer's datasheet
/// cap.val = 2.3e-6;
/// // Prints "Capacitor { val: 2300.000nF, esr: 40.0mΩ, esl: 0.4nH, lmnt: 1.0nH, package: SMD_2012M }"
/// println!("{cap}");
/// ```
#[derive(Debug, Copy, Clone, PartialEq)]
pub struct Capacitor {
    /// Capacitance in Farads
    pub val: f64,
    /// Equivalent Series Resistance
    pub esr: f64,
    /// Equivalent Series Inductance
    pub esl: f64,
    /// Mounting inductance
    pub lmnt: f64,
    /// Package type (e.g. 1005M or Eletrolytic)
    pub package: Package,
}

impl Capacitor {
    /// Calculate the resonant frequency formed by the capacitance in ESL+Lmnt
    /// 
    /// # Examples
    /// ```
    /// use complextest::passives::capacitor::Capacitor;
    /// use complextest::passives::packages::Package;
    /// 
    /// let cap = Capacitor{
    ///     val:3.0e-6,   //F
    ///     esr:0.050,    //Ω
    ///     esl:2.0e-9,   //H
    ///     lmnt:17.0e-9, //H
    ///     package:Package::SMD_1005M
    /// };
    /// let resonant = cap.resonant();
    /// let expected = 666626.0; //Hz
    /// let err = 1.0; //Hz
    /// println!("Resonant = {resonant:.1}Hz | expected = {expected:.1}Hz");
    /// assert!(resonant < (expected+err) && resonant > (expected-err));
    /// ```
    pub fn resonant(&self) -> f64 {
        1.0/(self.val*(self.esl+self.lmnt)).sqrt()/2.0/PI
    }

    /// Generate the equivalent capacitor impedance model in the s-domain
    /// 
    /// # Examples
    /// ```
    /// use complextest::passives::capacitor::Capacitor;
    /// use complextest::passives::packages::Package;
    /// use complextest::sdomain::Fs;
    /// 
    /// let cap = Capacitor{
    ///     val:1.0e-6,  //F
    ///     esr:0.089,   //Ω
    ///     esl:0.4e-9,  //H
    ///     lmnt:1.0e-9, //H
    ///     package:Package::Electrolytic
    /// };
    /// assert_eq!(
    ///     "( 1.000e0 + 8.900e-8s + 1.400e-15s² ) / ( 0.000e0 + 1.000e-6s )",
    ///     format!("{}", cap.model())
    /// )
    /// ```
    pub fn model(&self) -> Fs {
        sdomain::gen::rcl(self.esr, self.val, self.esl+self.lmnt)
    }

    /// Generate a capacitor model by searching for a template with matching value and package type.
    /// 
    /// # Arguments
    /// * `val` - Capacitance value in Farads to match
    /// * `package` - String representing the pacakge type to look for. For example: "0402" or "1005M"
    /// 
    /// # Examples
    /// ```
    /// use complextest::passives::capacitor::Capacitor;
    /// 
    /// let cap1 = Capacitor::from(4.7e-6, "0805");
    /// assert_eq!(
    ///     "Capacitor { val: 4700.000nF, esr: 4.0mΩ, esl: 0.5nH, lmnt: 1.0nH, package: SMD_2012M }",
    ///     format!("{}", cap1)
    /// );
    /// let cap2 = Capacitor::from(4.7e-6, "2012M");
    /// assert_eq!(
    ///     "Capacitor { val: 4700.000nF, esr: 4.0mΩ, esl: 0.5nH, lmnt: 1.0nH, package: SMD_2012M }",
    ///     format!("{}", cap2)
    /// );
    /// ```
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

    /// Generate a capacitor model by search for a template whose resonant frequency is within error
    /// of a given target.
    /// Returns an `Option<Capacitor>`, where None means no capacitor template was found whose resonant
    /// frequency was within error of the target.
    /// 
    /// # Arguments
    /// * `resonant` - The resonant frequency in Hertz to find the closest of
    /// * `err` - Acceptable error in Hertz with respect to the target resonant frequency
    /// 
    /// # Examples
    /// ```
    /// use complextest::passives::capacitor::Capacitor;
    /// 
    /// let cap = Capacitor::from_resonant(28.677e6, 0.01e6).unwrap();
    /// assert_eq!(
    ///     "Capacitor { val: 22.000nF, esr: 43.0mΩ, esl: 0.4nH, lmnt: 1.0nH, package: SMD_1005M }",
    ///     format!("{}", cap)
    /// );
    /// ```
    pub fn from_resonant(resonant: f64, err: f64) -> Option<Capacitor> {
        let database = include!("capacitor_database.in");
        let result = database
            .iter()
            .filter(|c| {
                let test_resonant = c.resonant();
                test_resonant < (resonant+err) && test_resonant > (resonant-err)
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
        assert_eq!("( 1.000e0 + 2.000e0s + 7.000e0s² ) / ( 0.000e0 + 1.000e0s )", format!("{}", cap.model()))
    }

    #[test]
    fn from_test() {
        let cap1 = Capacitor::from(4.7e-6, "0805");
        assert_eq!(
            "Capacitor { val: 4700.000nF, esr: 4.0mΩ, esl: 0.5nH, lmnt: 1.0nH, package: SMD_2012M }",
            format!("{}", cap1)
        );
        let cap2 = Capacitor::from(4.7e-6, "2012M");
        assert_eq!(
            "Capacitor { val: 4700.000nF, esr: 4.0mΩ, esl: 0.5nH, lmnt: 1.0nH, package: SMD_2012M }",
            format!("{}", cap2)
        );
    }

    #[test]
    fn from_resonant_test() {
        let cap1 = Capacitor::from_resonant(28.677e6, 0.01e6).unwrap();
        assert_eq!(
            "Capacitor { val: 22.000nF, esr: 43.0mΩ, esl: 0.4nH, lmnt: 1.0nH, package: SMD_1005M }",
            format!("{}", cap1)
        );
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