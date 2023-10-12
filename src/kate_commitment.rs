use crate::inner_product::Polynomial;

use zkstd::behave::{CurveGroup, FftField, Group, Pairing};

pub(crate) struct Proof<P: Pairing> {
    a: P::G2Affine,
    b: P::G1Affine,
    c: P::G1Affine,
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

    pub(crate) fn verify(&self, proof: Proof<P>) -> bool {
        let Proof { a, b, c } = proof;
        let lhs = P::pairing(b, a);
        let rhs = P::pairing(c, self.h);
        lhs == rhs
    }
}

#[cfg(test)]
mod tests {
    use super::{KateCommitment, Polynomial, Proof};

    use bls_12_381::{Fr as Scalar, G1Projective as G1, G2Projective as G2};
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

        // verifier sampling
        let b = Scalar::random(OsRng);

        // prover evaluation
        // evaluate with f(b)
        let c = poly.evaluate(b);
        // compute quotient polynomial f(x) - f(b) / x - c
        let q_poly = poly.divide(&c);
        // commit quotient polynomial
        let j = pp.commit(&q_poly);
        // x - c
        let g2_c = c * G2::ADDITIVE_GENERATOR;
        let h = pp.get_h();
        let i = (h - g2_c).into();
        // f(x) - c
        let p_c = pp.commit(&poly);
        let g1_c = c * G1::ADDITIVE_GENERATOR;
        let k = (p_c - g1_c).into();

        let proof = Proof { a: i, b: j, c: k };

        assert!(pp.verify(proof))
    }
}
