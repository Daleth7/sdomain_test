pub mod range_generators;
pub mod complex;
pub mod sdomain;
pub mod passives;
pub mod pdn;

#[cfg(test)]
mod tests {
    use super::*;

    use complex::Complex;
    use sdomain::{Fs, self};
    use rand::Rng;
    use std::cmp::max;


    fn print_sdomain(fs: &Fs) {
        let numer_str = fs.to_str_numer(false);
        let denom_str = fs.to_str_denom(false);
        let width = max(numer_str.len(), denom_str.len());
        println!("{numer_str:^width$}");
        println!("{:-^width$}", "");
        println!("{denom_str:^width$}");
    }

    #[test]
    fn miscellaneous_tests() {
        // Should not panic
        
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
        print_sdomain(&sdomain::gen::s_prec(1));
        println!("");
        print_sdomain(&sdomain::gen::step_prec(1));
        println!("");
        print_sdomain(&sdomain::gen::resistor_prec(10e3, 1));
        println!("");
        print_sdomain(&sdomain::gen::capacitor_prec(22e-6, 1));
        println!("");
        print_sdomain(&sdomain::gen::inductor_prec(0.47e-6, 1));
        println!("");
        print_sdomain(&sdomain::gen::rcl_prec(10e3, 22e-6, 0.47e-6, 1));
    }
}