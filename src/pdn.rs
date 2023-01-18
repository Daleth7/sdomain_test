use std::collections::HashMap;

use crate::sdomain::{Fs, parallel, parallel_multi};

/// Represents a Power Distribution Network (PDN) after a buck converter made up of the inductor and multiple output capacitors in parallel.
/// ```text
/// ┌─┤Inductor├────┬───────────┬─── ... ───┬─── Vout
/// │             ──┴──       ──┴──       ──┴──
/// │             Cap 0       Cap 1       Cap N
/// │             ──┬──       ──┬──       ──┬──
/// └───────────────┴───────────┴─── ... ───┴──┐
///                                            ⏚
/// ```
pub struct PDNModel {
    /// Inductor impedance model
    pub inductor: Fs,
    /// A collection of capacitor impedance models. Stored in a map where the key is a user-defined string.
    /// In addition to the model, the map also stores the quantity of each capacitor the PDN has.
    pub capacitors: HashMap<String, (u32, Fs)>,
}


impl PDNModel {
    /// Generate a new PDN from an inductor and a collection of capacitors
    /// 
    /// # Arguments
    /// * `inductor` - Inductor impedance model
    /// * `capacitors` - (Optional) A collection of capacitor impedance models. Expected
    ///                             to be a key value pair where:
    ///   * key = String
    ///   * value = (quantity, capacitor impedance model)
    /// 
    /// # Examples
    /// ```
    /// use sdomain_test::pdn::PDNModel;
    /// use sdomain_test::sdomain::gen;
    /// use sdomain_test::passives::capacitor::Capacitor;
    /// 
    /// // Create a PDN with just the inductor and no output caps. Useful if the
    /// // output caps. will be added later.
    /// let mut pdn = PDNModel::from(gen::rl(0.025/*Ω*/, 330e-9/*H*/), None);
    /// 
    /// // Add two of the same output capacitor to the network
    /// pdn.add_capacitor("89μF (derated) ceramic (2012M)", gen::rcl(0.003/*Ω*/, 89e-6/*F*/, 1.2e-9/*H*/), 2);
    /// ```
    pub fn from(inductor: Fs, capacitors: Option<HashMap<String, (u32, Fs)>>) -> Self {
        let cap_collection = match capacitors {
            Some(collection) => collection,
            None => HashMap::<String, (u32, Fs)>::new(),
        };
        Self{inductor, capacitors: cap_collection}
    }

    /// Add a new output capacitor to the network.
    /// 
    /// # Arguments
    /// * `id` - Identifier to give to capacitor. Useful if later on the capacitor needs to be
    ///        extracted from the network or if the quantity needs to be changed.
    ///        If the network already contains a capacitor with this id, the capacitor's
    ///        quantity is updated.
    /// * `capacitor` - Capacitor impedance model to add.
    /// * `qty` - Quantity.
    /// 
    /// # Examples
    /// ```
    /// use sdomain_test::pdn::PDNModel;
    /// use sdomain_test::sdomain::gen;
    /// use sdomain_test::passives::capacitor::Capacitor;
    /// 
    /// let mut pdn = PDNModel::from(gen::rl(0.025/*Ω*/, 330e-9/*H*/), None);
    /// 
    /// // Add two of the same output capacitor to the network
    /// pdn.add_capacitor("89μF (derated) ceramic (2012M)", gen::rcl(0.003/*Ω*/, 89e-6/*F*/, 1.2e-9/*H*/), 2);
    /// 
    /// // Add capacitor generated from template
    /// let hf_cap = Capacitor::from(22e-9, "0201"); // Or "0603M" for metric
    /// pdn.add_capacitor("High Frequency Cap", hf_cap.model(), 1);
    pub fn add_capacitor(&mut self, id: &str, capacitor: Fs, qty: u32) {
        let new_key = id.to_string();
        if self.capacitors.contains_key(&new_key) {
            self.capacitors.get_mut(&new_key).unwrap().0 += qty;
        } else {
            self.capacitors.insert(new_key, (qty, capacitor));
        }
    }

    /// Generate an impedance model of the PDN in the s-domain.
    /// 
    /// # Examples
    /// ```
    /// use sdomain_test::pdn::PDNModel;
    /// use sdomain_test::sdomain::gen;
    /// use sdomain_test::passives::capacitor::Capacitor;
    /// 
    /// let mut pdn = PDNModel::from(gen::rl(0.025/*Ω*/, 330e-9/*H*/), None);
    /// pdn.add_capacitor("89μF (derated) ceramic (2012M)", gen::rcl(0.003/*Ω*/, 89e-6/*F*/, 1.2e-9/*H*/), 2);
    /// let hf_cap = Capacitor::from(22e-9, "0201"); // Capacitor { val: 22.000nF, esr: 39.0mΩ, esl: 0.3nH, lmnt: 1.0nH, package: SMD_0603M }
    /// pdn.add_capacitor("High Frequency Cap", hf_cap.model(), 1);
    /// 
    /// // Generate impedance model
    /// let Zpdn = pdn.model();
    /// assert_eq!(
    ///     "( 1.250e-2 + 1.683e-7s + 4.553e-14s² + 1.767e-20s³ + 1.642e-29s⁴ + 5.040e-37s⁵ ) / ( 5.000e-1 + 2.359e-6s + 2.943e-11s² + 2.631e-20s³ + 1.229e-27s⁴ )",
    ///     format!("{}", Zpdn)
    /// );
    pub fn model(&self) -> Fs {
        // Collapse collection of components into one equation
        let mul_caps = self.capacitors.iter().map(|(_, (qty, fs))| fs.clone().scale(1.0 / *qty as f64));
        parallel(parallel_multi(mul_caps), self.inductor.clone())
        // let zrm = self.inductor.clone();
        // match parallel_multi(mul_caps) {
        //     Some(zout) => parallel(zout, zrm),
        //     None => zrm,
        // }
    }
}