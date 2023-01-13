use std::collections::HashMap;

use crate::sdomain::{Fs, self};

pub struct PDNModel {
    pub inductor: Fs,
    pub capacitors: HashMap<String, (u32, Fs)>,
}


impl PDNModel {
    pub fn new() -> Self {
        Self{inductor: sdomain::gen::inductor(0.47e-6), capacitors: HashMap::<String, (u32, Fs)>::new()}
    }

    pub fn from(inductor: Fs, capacitors: Option<HashMap<String, (u32, Fs)>>) -> Self {
        let cap_collection = match capacitors {
            Some(collection) => collection,
            None => HashMap::<String, (u32, Fs)>::new(),
        };
        Self{inductor, capacitors: cap_collection}
    }

    pub fn add_capacitor(&mut self, id: String, capacitor: Fs, qty: u32) {
        if self.capacitors.contains_key(&id) {
            self.capacitors.get_mut(&id).unwrap().0 += qty;
        } else {
            self.capacitors.insert(id, (qty, capacitor));
        }
    }

    pub fn model(&self) -> Fs {
        // Collapse collection of components into one equation
        let mul_caps: Vec<Fs> = self.capacitors.iter().map(|(_, (qty, fs))| fs.clone().scale(1.0 / *qty as f64)).collect();
        let mut cap_iter = mul_caps.into_iter();
        match cap_iter.next() {
            Some(init_cap) => {
                let mut go = init_cap.invert();
                for cap in cap_iter {
                    go += &cap.invert();
                }
                (go + &self.inductor.clone().invert()).invert()
            }
            None => self.inductor.clone()
        }
    }
}