use crate::inner_product::Polynomial;

use zkstd::behave::{CurveGroup, FftField, Group, Pairing};

pub(crate) struct KateCommitment<P: Pairing> {
    g: Vec<P::G1Affine>,
    h: P::G2Affine,
}

impl<P: Pairing> KateCommitment<P> {
    pub(crate) fn new(k: usize, r: P::ScalarField) -> Self {
        // G1, r * G1, r^2 * G1, ..., r^n-1 * G1
        let g = (0..=1 << k)
            .map(|i| {
                let tw = P::G1Projective::ADDITIVE_GENERATOR * r.pow(i);
                P::G1Affine::from(tw)
            })
            .collect::<Vec<_>>();
        let h = (P::G2Affine::ADDITIVE_GENERATOR * r).into();
        Self { g, h }
    }

    /// c_0 + c_1 * x + c_2 * x^2 + ... + c_d * x^d
    pub(crate) fn commit(&self, polynomial: &Polynomial<P::ScalarField>) -> P::G1Affine {
        polynomial
            .coeffs
            .iter()
            .zip(self.g.iter())
            .fold(P::G1Projective::ADDITIVE_IDENTITY, |sum, (scalar, base)| {
                sum + *base * *scalar
            })
            .into()
    }
}

#[cfg(test)]
mod tests {
    use super::{KateCommitment, Polynomial};

    use bls_12_381::{Fr as Scalar, G1Projective as Point};
    use ec_pairing::TatePairing;
    use rand::rngs::OsRng;
    use zkstd::{
        behave::{CurveGroup, Group},
        common::CurveAffine,
    };

    #[test]
    fn commit_test() {
        let k = 8;
        let r = Scalar::random(OsRng);
        let pp = KateCommitment::<TatePairing>::new(k, r);
        let coeffs = (0..1 << k).map(|_| Scalar::random(OsRng)).collect();
        let polynomial = Polynomial::new(coeffs);
        let commitment = pp.commit(&polynomial);
        let eval = polynomial.evaluate(r);

        assert_eq!(commitment.to_extended(), eval * Point::ADDITIVE_GENERATOR)
    }
}
