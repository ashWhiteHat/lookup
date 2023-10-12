use crate::inner_product::Polynomial;

use rayon::join;
use zkstd::common::FftField;

pub struct Fft<F: FftField> {
    // domain size
    n: usize,
    // n th root of unity
    twiddle_factors: Vec<F>,
    // n th root of unity inverse
    inv_twiddle_factors: Vec<F>,
    // n inverse
    n_inv: F,
    // bit reverse index
    bit_reverse: Vec<(usize, usize)>,
}

impl<F: FftField> Fft<F> {
    pub fn new(k: usize) -> Self {
        assert!(k >= 1);
        let n = 1 << k;
        let half_n = n / 2;
        let offset = 64 - k;

        // compute twiddle factors
        let g = (0..F::S - k).fold(F::ROOT_OF_UNITY, |acc, _| acc.square());
        let twiddle_factors = (0..half_n)
            .scan(F::one(), |w, _| {
                let tw = *w;
                *w *= g;
                Some(tw)
            })
            .collect::<Vec<_>>();

        // compute inverse twiddle factors
        let g_inv = g.invert().unwrap();
        let inv_twiddle_factors = (0..half_n)
            .scan(F::one(), |w, _| {
                let tw = *w;
                *w *= g_inv;
                Some(tw)
            })
            .collect::<Vec<_>>();

        let bit_reverse = (0..n as u64)
            .filter_map(|i| {
                let r = i.reverse_bits() >> offset;
                (i < r).then_some((i as usize, r as usize))
            })
            .collect::<Vec<_>>();

        let n_inv = F::from(n as u64).invert().unwrap();

        Self {
            n,
            twiddle_factors,
            inv_twiddle_factors,
            n_inv,
            bit_reverse,
        }
    }

    /// perform discrete fourier transform
    pub(crate) fn dft(&self, poly: &mut Polynomial<F>) {
        self.prepare_fft(poly);
        classic_fft_arithmetic(&mut poly.coeffs, self.n, 1, &self.twiddle_factors)
    }

    /// perform classic inverse discrete fourier transform
    pub(crate) fn idft(&self, poly: &mut Polynomial<F>) {
        self.prepare_fft(poly);
        classic_fft_arithmetic(&mut poly.coeffs, self.n, 1, &self.inv_twiddle_factors);
        poly.coeffs
            .iter_mut()
            .for_each(|coeff| *coeff *= self.n_inv)
    }

    /// polynomial multiplication
    pub(crate) fn poly_mul(&self, mut rhs: Polynomial<F>, mut lhs: Polynomial<F>) -> Polynomial<F> {
        self.dft(&mut rhs);
        self.dft(&mut lhs);
        let mut mul_poly = Polynomial::new(
            rhs.coeffs
                .iter()
                .zip(lhs.coeffs.iter())
                .map(|(a, b)| *a * *b)
                .collect(),
        );
        self.idft(&mut mul_poly);
        mul_poly
    }

    fn prepare_fft(&self, poly: &mut Polynomial<F>) {
        poly.coeffs.resize(self.n, F::zero());
        self.bit_reverse
            .iter()
            .for_each(|(i, ri)| poly.coeffs.swap(*ri, *i));
    }
}

fn classic_fft_arithmetic<F: FftField>(
    coeffs: &mut [F],
    n: usize,
    twiddle_chunk: usize,
    twiddles: &[F],
) {
    if n == 2 {
        let t = coeffs[1];
        coeffs[1] = coeffs[0];
        coeffs[0] += t;
        coeffs[1] -= t;
    } else {
        let (left, right) = coeffs.split_at_mut(n / 2);
        join(
            || classic_fft_arithmetic(left, n / 2, twiddle_chunk * 2, twiddles),
            || classic_fft_arithmetic(right, n / 2, twiddle_chunk * 2, twiddles),
        );
        butterfly_arithmetic(left, right, twiddle_chunk, twiddles)
    }
}

fn butterfly_arithmetic<F: FftField>(
    left: &mut [F],
    right: &mut [F],
    twiddle_chunk: usize,
    twiddles: &[F],
) {
    // case when twiddle factor is one
    let t = right[0];
    right[0] = left[0];
    left[0] += t;
    right[0] -= t;

    left.iter_mut()
        .zip(right.iter_mut())
        .enumerate()
        .skip(1)
        .for_each(|(i, (a, b))| {
            let mut t = *b;
            t *= twiddles[i * twiddle_chunk];
            *b = *a;
            *a += t;
            *b -= t;
        });
}

#[cfg(test)]
mod tests {
    use super::*;

    use bls_12_381::Fr as Scalar;
    use rand::rngs::OsRng;
    use zkstd::behave::{Group, PrimeField};

    fn arb_poly(k: u32) -> Polynomial<Scalar> {
        Polynomial {
            coeffs: (0..(1 << k))
                .map(|_| Scalar::random(OsRng))
                .collect::<Vec<_>>(),
        }
    }

    fn naive_multiply<F: PrimeField>(a: Vec<F>, b: Vec<F>) -> Vec<F> {
        assert_eq!(a.len(), b.len());
        let mut c = vec![F::zero(); a.len() + b.len()];
        a.iter().enumerate().for_each(|(i_a, coeff_a)| {
            b.iter().enumerate().for_each(|(i_b, coeff_b)| {
                c[i_a + i_b] += *coeff_a * *coeff_b;
            })
        });
        c
    }

    fn point_mutiply<F: FftField>(a: Polynomial<F>, b: Polynomial<F>) -> Polynomial<F> {
        assert_eq!(a.coeffs.len(), b.coeffs.len());
        Polynomial {
            coeffs: a
                .coeffs
                .iter()
                .zip(b.coeffs.iter())
                .map(|(coeff_a, coeff_b)| *coeff_a * *coeff_b)
                .collect::<Vec<F>>(),
        }
    }

    #[test]
    fn fft_transformation_test() {
        let mut poly_a = arb_poly(10);
        let poly_b = poly_a.clone();
        let classic_fft = Fft::new(10);

        classic_fft.dft(&mut poly_a);
        classic_fft.idft(&mut poly_a);

        assert_eq!(poly_a, poly_b)
    }

    #[test]
    fn fft_multiplication_test() {
        let coeffs_a = arb_poly(4);
        let coeffs_b = arb_poly(4);
        let fft = Fft::new(5);
        let poly_c = coeffs_a.clone();
        let poly_d = coeffs_b.clone();
        let mut poly_a = coeffs_a;
        let mut poly_b = coeffs_b;
        let poly_g = poly_a.clone();
        let poly_h = poly_b.clone();

        let poly_e = Polynomial {
            coeffs: naive_multiply(poly_c.coeffs, poly_d.coeffs),
        };

        fft.dft(&mut poly_a);
        fft.dft(&mut poly_b);
        let mut poly_f = point_mutiply(poly_a, poly_b);
        fft.idft(&mut poly_f);

        let poly_i = fft.poly_mul(poly_g, poly_h);

        assert_eq!(poly_e, poly_f);
        assert_eq!(poly_e, poly_i)
    }

    #[test]
    fn inverse_fft_multiplication_test() {
        let coeffs_a = arb_poly(4);
        let coeffs_b = arb_poly(4);
        let fft = Fft::new(5);
        let poly_c = coeffs_a.clone();
        let poly_d = coeffs_b.clone();
        let mut poly_a = coeffs_a;
        let mut poly_b = coeffs_b;
        let poly_g = poly_a.clone();
        let poly_h = poly_b.clone();

        let poly_e = Polynomial {
            coeffs: naive_multiply(poly_c.coeffs, poly_d.coeffs),
        };

        fft.dft(&mut poly_a);
        fft.dft(&mut poly_b);
        let mut poly_f = point_mutiply(poly_a, poly_b);
        fft.idft(&mut poly_f);

        let poly_i = fft.poly_mul(poly_g, poly_h);

        assert_eq!(poly_e, poly_f);
        assert_eq!(poly_e, poly_i)
    }
}
