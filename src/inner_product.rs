use rand::rngs::OsRng;
use zkstd::common::FftField;

use crate::fft::Fft;

// first to last, x^0 to x^n-1
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct Polynomial<F: FftField> {
    pub(crate) coeffs: Vec<F>,
}

impl<F: FftField> Polynomial<F> {
    pub(crate) fn new(coeffs: Vec<F>) -> Self {
        Self { coeffs }
    }

    // log(n) inner product
    pub(crate) fn inner_product(&self, rhs: &Self) -> F {
        assert_eq!(self.coeffs.len(), rhs.coeffs.len());
        self.coeffs
            .iter()
            .zip(rhs.coeffs.iter())
            .fold(F::zero(), |sum, (a, b)| sum + *a * *b)
    }

    pub(crate) fn random(k: usize) -> Self {
        let n = 1 << k;
        let coeffs = (0..n).map(|_| F::random(OsRng)).collect();
        Self { coeffs }
    }

    pub(crate) fn half(self) -> (Self, Self) {
        let n = self.coeffs.len();
        let half = n / 2;
        let (lo, hi) = self.coeffs.split_at(half);
        (
            Self {
                coeffs: lo.to_vec(),
            },
            Self {
                coeffs: hi.to_vec(),
            },
        )
    }

    pub(crate) fn scalar(self, scalar: F) -> Self {
        let coeffs = self.coeffs.iter().map(|coeff| *coeff * scalar).collect();
        Self { coeffs }
    }

    pub(crate) fn evaluate(&self, at: F) -> F {
        let mut acc = F::one();
        self.coeffs.iter().fold(F::zero(), |sum, coeff| {
            let tmp = acc;
            acc *= at;
            sum + *coeff * tmp
        })
    }

    pub(crate) fn divide(&self, at: &F) -> Self {
        let mut coeffs = self
            .coeffs
            .iter()
            .rev()
            .scan(F::zero(), |w, coeff| {
                let tmp = *w + coeff;
                *w = tmp * at;
                Some(tmp)
            })
            .collect::<Vec<_>>();
        coeffs.pop();
        coeffs.reverse();
        Self { coeffs }
    }

    pub(crate) fn from_evals(evals: Vec<F>, fft: Fft<F>) -> Self {
        let mut coeffs = Self { coeffs: evals };
        fft.idft(&mut coeffs);
        coeffs
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    use bls_12_381::Fr as Scalar;
    use zkstd::common::{Group, PrimeField};

    fn poly_mul<F: FftField>(a: &Polynomial<F>, b: &Polynomial<F>) -> Polynomial<F> {
        let size = a.coeffs.len() + b.coeffs.len() - 1;
        let mut coeffs = (0..size).map(|_| F::zero()).collect::<Vec<_>>();
        for (i, ac) in a.coeffs.iter().enumerate() {
            for (j, bc) in b.coeffs.iter().enumerate() {
                coeffs[i + j] += *ac * *bc;
            }
        }
        Polynomial { coeffs }
    }

    #[test]
    fn inner_product_proof_test() {
        // setup
        let k = 8;
        let x = Scalar::random(OsRng);
        let x_inv = x.invert().unwrap();
        let xx = x.square();
        let xx_inv = x_inv.square();
        let a_poly = Polynomial::<Scalar>::random(k);
        let b_poly = Polynomial::<Scalar>::random(k);

        // compress
        let (alo, ahi) = a_poly.half();
        let (blo, bhi) = b_poly.half();
        let c = alo.inner_product(&blo) + ahi.inner_product(&bhi);
        let l = alo.inner_product(&bhi);
        let r = ahi.inner_product(&blo);
        let xx_l = xx * l;
        let xx_inv_r = xx_inv * r;
        let c_prime = c + xx_l + xx_inv_r;
    }

    #[test]
    fn inner_product_test() {
        let k = 8;
        let a_poly = Polynomial::<Scalar>::random(k);
        let b_poly = Polynomial::<Scalar>::random(k);
        let naive_product = a_poly.inner_product(&b_poly);

        let (alo, ahi) = a_poly.half();
        let (blo, bhi) = b_poly.half();
        let half_product = alo.inner_product(&blo) + ahi.inner_product(&bhi);

        assert_eq!(naive_product, half_product)
    }

    #[test]
    fn multiplication_test() {
        let k = 3;
        let r = Scalar::random(OsRng);
        let coeffs = (0..1 << k)
            .map(|_| Scalar::random(OsRng))
            .collect::<Vec<_>>();
        let factor_coeffs = vec![-r, Scalar::one()];
        let a = Polynomial { coeffs };
        let b = Polynomial {
            coeffs: factor_coeffs,
        };
        // a * b
        let c = poly_mul(&a, &b);
        // a * b / b
        let d = c.divide(&r);

        assert_eq!(a, d);
    }
}
