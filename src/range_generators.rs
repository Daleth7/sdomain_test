use std::ops::{Mul, MulAssign, Div};

/// Generate a range of values logramithically.
/// 
/// # Arguments
/// * `start`    - Starting value in the range
/// * `end`      - End of the range, non-inclusive
/// * `radix`    - Number base
/// * `division` - Amount of numbers to generate for each power of radix
/// 
/// # Examples
/// ```
/// use complextest::range_generators::gen_log_range;
/// 
/// let log10_range = gen_log_range(1, 1000, 10, 10);
/// assert_eq!(
///     vec![
///         1, 2, 3, 4, 5, 6, 7, 8, 9,
///         10, 20, 30, 40, 50, 60, 70, 80, 90,
///         100, 200, 300, 400, 500, 600, 700, 800, 900
///     ],
///     log10_range
/// );
/// 
/// // Use a different number base and generate floating points instead of integers
/// let log3f_range = gen_log_range(1.0, 30.0, 3.0, 3);
/// assert_eq!("[ 1.0, 2.0, 3.0, 6.0, 9.0, 18.0, 27.0, ]", vec2str(&log3f_range));
/// 
/// // Use more divisions per power
/// let log10f_range = gen_log_range(1.0, 25.0, 10.0, 20);
/// assert_eq!("[ 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0, 15.0, 20.0, ]", vec2str(&log10f_range));
/// 
/// fn vec2str(collection: &Vec<f64>) -> String {
///     let mut result = "[ ".to_string();
///     for n in collection {
///         result.push_str(&format!("{:.1}, ", n));
///     }
///     result.push_str("]");
///     result
/// }
/// ```
pub fn gen_log_range<T>(mut start: T, end: T, radix: T, divisions: i32) -> Vec<T> where
    T: Mul<Output = <T as Div>::Output> + MulAssign + Div + PartialOrd + From<i32> + Copy,
    <T as Div>::Output: Mul,
    <T as Div>::Output: Copy,
    Vec<T>: FromIterator<<<T as Mul>::Output as Mul<<T as Div>::Output>>::Output>,
    <<T as Mul>::Output as Mul<<T as Div>::Output>>::Output: PartialOrd<T>
{
    let gain = radix / T::from(divisions);
    let mut toreturn = Vec::<T>::new();
    while start < end {
        toreturn.extend_from_slice(
            &Vec::<T>::from((1..divisions)
                .map(|x| T::from(x)*start*gain)
                .filter(|x| *x < end && *x >= start)
                .collect::<Vec<T>>())
        );
        start *= radix.clone();
    }
    toreturn
}


#[cfg(test)]
mod tests {
    use super::*;

    fn gen_exp_log10() -> Vec<i32> {
        vec![
            1  , 2  , 3  , 4  , 5  , 6  , 7  , 8  , 9  ,
            10 , 20 , 30 , 40 , 50 , 60 , 70 , 80 , 90 ,
            100, 200, 300, 400, 500, 600, 700, 800, 900,
            ]
    }

    fn gen_exp_log7() -> Vec<i32> {
        vec![
            1  , 2  , 3  , 4  , 5  , 6  ,
            7  , 14 , 21 , 28 , 35 , 42 ,
            49 , 98 , 147, 196
            ]
    }

    #[test]
    fn gen_log_range_i32_test() {
        assert_eq!(gen_exp_log10(), gen_log_range(1, 1000, 10, 10));
        assert_eq!(gen_exp_log7(), gen_log_range(1, 200, 7, 7));
    }

    #[test]
    fn gen_log_range_f64_test() {
        let exp_range_radix10 = gen_exp_log10()
            .iter()
            .map(|x| *x as f64)
            .collect::<Vec<f64>>();
        let exp_range_radix7 = gen_exp_log7()
            .iter()
            .map(|x| *x as f64)
            .collect::<Vec<f64>>();
        assert_eq!(exp_range_radix10, gen_log_range(1.0, 1000.0, 10.0, 10));
        assert_eq!(exp_range_radix7, gen_log_range(1.0, 200.0, 7.0, 7));
    }

    #[test]
    fn gen_log_range_f64_stepsize_test() {
        let exp_range_radix10 = vec![
            1.0 , 1.5 , 2.0 , 2.5 , 3.0 , 3.5 , 4.0 , 4.5 , 5.0 , 5.5 , 6.0 , 6.5 , 7.0 , 7.5 , 8.0 , 8.5 , 9.0 , 9.5,
            10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 50.0, 55.0, 60.0
            ];
        let exp_range_radix4 = vec![
            1.00, 1.25, 1.50, 1.75, 2.00, 2.25, 2.50 , 2.75 , 3.00 , 3.25 , 3.50 , 3.75,
            4.00, 5.00, 6.00, 7.00, 8.00, 9.00, 10.00, 11.00, 12.00, 13.00, 14.00
            ];
        assert_eq!(exp_range_radix10, gen_log_range(1.0, 61.0, 10.0, 20));
        assert_eq!(exp_range_radix4, gen_log_range(1.0, 14.5, 4.0, 16));
    }
}