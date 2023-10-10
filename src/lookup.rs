use zkstd::common::PrimeField;

pub(crate) struct XORTable<F: PrimeField> {
    n: usize,
    a: Vec<F>,
    b: Vec<F>,
    c: Vec<F>,
}

impl<F: PrimeField> XORTable<F> {
    pub(crate) fn precompute() -> Self {
        let bit_length = 4;
        let n = 1 << bit_length;
        let (mut a, mut b, mut c) = (Vec::new(), Vec::new(), Vec::new());
        for i in 0..n {
            for j in 0..n {
                let k = i ^ j;
                a.push(F::from(i as u64));
                b.push(F::from(j as u64));
                c.push(F::from(k as u64));
            }
        }
        Self { n, a, b, c }
    }

    pub(crate) fn compress(&self, α: F) -> Vec<F> {
        let αα = α.square();
        self.a
            .iter()
            .zip(self.b.iter())
            .zip(self.c.iter())
            .map(|((t1, t2), t3)| *t1 + α * t2 + αα * t3)
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::XORTable;
    use bls_12_381::Fr as Scalar;
    use rand::rngs::OsRng;
    use rand::{thread_rng, Rng};
    use zkstd::common::{Group, PrimeField};

    fn compress_witness<F: PrimeField>(a: u64, b: u64, c: u64, α: F) -> F {
        F::from(a) + α * F::from(b) + α.square() * F::from(c)
    }

    #[test]
    fn table_generation_test() {
        let bit_length = 4;
        let length = 1 << (bit_length * 2);
        let xor_table = XORTable::<Scalar>::precompute();
        assert_eq!(xor_table.a.len(), length);
        assert_eq!(xor_table.b.len(), length);
        assert_eq!(xor_table.c.len(), length);
        for ((i, j), k) in xor_table
            .a
            .iter()
            .zip(xor_table.b.iter())
            .zip(xor_table.c.iter())
        {
            assert_eq!(i ^ j, *k)
        }
    }

    #[test]
    fn lookup_test() {
        let bit_length = 4;
        let range = 1 << bit_length;
        let a = thread_rng().gen_range(0..range);
        let b = thread_rng().gen_range(0..range);
        let c = a ^ b;
        let d = a + b;
        let α = Scalar::random(OsRng);
        let xor_table = XORTable::<Scalar>::precompute();
        let ti = xor_table.compress(α);
        let fi = compress_witness(a, b, c, α);
        let gi = compress_witness(a, b, d, α);

        assert!(ti.contains(&fi));
        assert!(!ti.contains(&gi))
    }
}
