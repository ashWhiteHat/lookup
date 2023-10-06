use zkstd::behave::{CurveAffine, Group};
use zkstd::common::RngCore;

pub(crate) struct PedersenCommitment<C: CurveAffine> {
    g: C,
    h: C,
}

impl<C: CurveAffine> PedersenCommitment<C> {
    pub(crate) fn new(r: impl RngCore) -> Self {
        let g = C::ADDITIVE_GENERATOR;
        let h = (g * C::Scalar::random(r)).into();
        Self { g, h }
    }

    pub(crate) fn commit(&self, s: C::Scalar, r: impl RngCore) -> (C, C::Scalar) {
        let t = C::Scalar::random(r);
        let commitment = self.g * s + self.h * t;
        (commitment.into(), t)
    }

    pub(crate) fn open(&self, commitment: C, s: C::Scalar, t: C::Scalar) -> bool {
        commitment.to_extended() == self.g * s + self.h * t
    }
}

#[cfg(test)]
mod tests {
    use super::PedersenCommitment;

    use bls_12_381::{Fr as Scalar, G1Affine as Point};
    use rand::rngs::OsRng;
    use zkstd::common::Group;

    #[test]
    fn perdersen_commitment_test() {
        let s = Scalar::random(OsRng);
        let params = PedersenCommitment::<Point>::new(OsRng);
        let (commitment, t) = params.commit(s, OsRng);
        assert!(params.open(commitment, s, t))
    }
}
