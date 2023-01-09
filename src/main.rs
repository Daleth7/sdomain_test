use crate::{complex::Complex, sdomain::Fs};
use rand::Rng;
use std::cmp::max;

mod complex;
mod sdomain;

fn print_sdomain(fs: &Fs) {
    let numer_str = fs.to_str_numer(false);
    let denom_str = fs.to_str_denom(false);
    let width = max(numer_str.len(), denom_str.len());
    println!("{numer_str:^width$}");
    println!("{:-^width$}", "");
    println!("{denom_str:^width$}");
}

fn main() {
    println!("Hello, world!");

    let mut prng = rand::thread_rng();

    for _ in 1..=100 {
        let c = Complex::new(prng.gen(), prng.gen());
        let mag = c.mag();
        let phase = c.phase();
        let phase_deg = c.phase_deg();
        println!("{c} (M: {mag:.3}, P: {phase:.3} ({phase_deg:>4.1}°))");
    }
    println!("\n\n");

    let fs = Fs{numer:vec![1.0, -2.3, 4.5], denom:vec![-0.67, 4.1], prec:3};
    print_sdomain(&fs);

    println!("");
    print_sdomain(&sdomain::gen::s(1));
    println!("");
    print_sdomain(&sdomain::gen::step(1));
    println!("");
    print_sdomain(&sdomain::gen::resistor(10e3, 1));
    println!("");
    print_sdomain(&sdomain::gen::capacitor(22e-6, 1));
    println!("");
    print_sdomain(&sdomain::gen::inductor(0.47e-6, 1));
    println!("");
    print_sdomain(&sdomain::gen::rcl(10e3, 22e-6, 0.47e-6, 1));
}
