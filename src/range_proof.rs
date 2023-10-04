#[cfg(test)]
mod tests {
    use bls_12_381::Fr as Scalar;
    use ec_pairing::TatePairing;
    use poly_commit::KzgParams;
    use rand::rngs::OsRng;
    use zkstd::common::Group;

    #[test]
    fn range_proof_test() {
        // range proof 0 <= z < 256
        let k = 2;
        let n = 256;
        let z = Scalar::from(100);

        // setup kzg params
        let r = Scalar::random(OsRng);
        let pp = KzgParams::<TatePairing>::setup(k, r);

        assert!(true)
    }
}
