use crate::inner_product::Polynomial;

use zkstd::behave::{CurveGroup, FftField, Pairing};

pub(crate) struct Proof<P: Pairing> {
    a: P::G1Affine,
    b: P::G2Affine,
    c: P::G1Affine,
}

impl<P: Pairing> Proof<P> {
    pub(crate) fn new(a: P::G1Affine, b: P::G2Affine, c: P::G1Affine) -> Self {
        Self { a, b, c }
    }

    pub(crate) fn verify(self) -> bool {
        let Proof { a, b, c } = self;
        let lhs = P::pairing(a, b);
        let rhs = P::pairing(c, P::G2Affine::ADDITIVE_GENERATOR);
        lhs == rhs
    }
}

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

    pub(crate) fn get_h(&self) -> P::G2Affine {
        self.h
    }
}

#[cfg(test)]
mod tests {
    use super::{KateCommitment, Polynomial, Proof};

    use bls_12_381::{Fr as Scalar, G1Affine as G1, G2Affine as G2};
    use ec_pairing::TatePairing;
    use rand::rngs::OsRng;
    use zkstd::behave::{CurveAffine, CurveGroup, Group, Pairing};

    fn sample_data<P: Pairing>(
        r: P::ScalarField,
    ) -> (Polynomial<P::ScalarField>, KateCommitment<P>) {
        let k = 8;
        let pp = KateCommitment::<P>::new(k, r);
        let coeffs = (0..1 << k).map(|_| P::ScalarField::random(OsRng)).collect();
        let poly = Polynomial::new(coeffs);
        (poly, pp)
    }

    #[test]
    fn commit_test() {
        let r = Scalar::random(OsRng);
        let (poly, pp) = sample_data::<TatePairing>(r);
        let commitment = pp.commit(&poly);
        let eval = poly.evaluate(r);

        assert_eq!(commitment.to_extended(), eval * G1::ADDITIVE_GENERATOR)
    }

    #[test]
    fn kzg_test() {
        // setup params
        let r = Scalar::random(OsRng);
        let (poly, pp) = sample_data::<TatePairing>(r);

        // verifier
        let b = Scalar::random(OsRng);

        // prover
        // 1. params computation

        // f(b): evaluate at b
        let b_eval = poly.evaluate(b);
        // q(x): compute quotient f(x) - f(b) / x - b
        let q_poly = poly.divide(&b);
        let b_g2 = G2::ADDITIVE_GENERATOR * b;
        let h = pp.get_h();

        // 2. generate proof

        // a - b
        let a = pp.commit(&q_poly);
        // commit q(a)
        let b = (h - b_g2).into();
        // f(a) - f(b)
        let c = (pp.commit(&poly) - G1::ADDITIVE_GENERATOR * b_eval).into();

        let proof: Proof<TatePairing> = Proof::new(a, b, c);

        // 3. proof verification
        assert!(proof.verify())
    }
}
