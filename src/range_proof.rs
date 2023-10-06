// range proof implementation describe as follows
// https://hackmd.io/@dabo/B1U4kx8XI#Range-proof-for-the-range-02n

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

// describe z as binary representation
fn z_to_g<F: PrimeField>(z: F) -> Vec<u8> {
    z.to_bits()
}

#[cfg(test)]
mod tests {
    use super::binary_check;

    use bls_12_381::Fr as Scalar;
    use rand::rngs::OsRng;
    use zkstd::common::Group;

    #[test]
    fn binary_check_test() {
        let z = Scalar::random(OsRng);
        assert!(binary_check(z))
    }
}
