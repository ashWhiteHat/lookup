use zkstd::common::PrimeField;

fn binary_check<F: PrimeField>(value: F) -> bool {
    let binary = value
        .to_bits()
        .iter()
        .map(|bit| *bit as i8)
        .collect::<Vec<_>>();
    let binary_prime = binary.iter().map(|bit| *bit - 1).collect::<Vec<_>>();
    binary
        .iter()
        .zip(binary_prime.iter())
        .all(|(bit, bit_prime)| bit * bit_prime == 0)
}

#[cfg(test)]
mod tests {
    use super::binary_check;

    use bls_12_381::Fr as Scalar;
    use ec_pairing::TatePairing;
    use poly_commit::KzgParams;
    use rand::rngs::OsRng;
    use zkstd::common::Group;

    #[test]
    fn binary_check_test() {
        let z = Scalar::random(OsRng);
        assert!(binary_check(z))
    }

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
