use crate::inner_product::Polynomial;

use zkstd::behave::{CurveGroup, FftField, Group, Pairing};
use zkstd::common::RngCore;

pub(crate) struct KateCommitment<P: Pairing> {
    g: Vec<P::G1Affine>,
    h: P::G2Affine,
}

impl<P: Pairing> KateCommitment<P> {
    pub(crate) fn new(k: usize, r: impl RngCore) -> Self {
        let t = P::ScalarField::random(r);
        let g = (0..=1 << k)
            .map(|i| {
                let tw = P::G1Projective::ADDITIVE_GENERATOR * t.pow(i);
                P::G1Affine::from(tw)
            })
            .collect::<Vec<_>>();
        let h = (P::G2Affine::ADDITIVE_GENERATOR * t).into();
        Self { g, h }
    }

    pub(crate) fn commit(&self, polynomial: Polynomial<P::ScalarField>) -> P::G1Affine {
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
