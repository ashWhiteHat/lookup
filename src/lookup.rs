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

    pub(crate) fn compress(&self, alpha: F) -> Vec<F> {
        let alpha2 = alpha.square();
        self.a
            .iter()
            .zip(self.b.iter())
            .zip(self.c.iter())
            .map(|((t1, t2), t3)| *t1 + alpha * t2 + alpha2 * t3)
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

    fn witness_vectors<F: PrimeField>(range: u64, alpha: F) -> Vec<F> {
        let i = 24;
        (0..i)
            .map(|_| {
                let a = thread_rng().gen_range(0..range);
                let b = thread_rng().gen_range(0..range);
                let c = a ^ b;
                compress_witness(a, b, c, alpha)
            })
            .collect()
    }

    fn compress_witness<F: PrimeField>(a: u64, b: u64, c: u64, alpha: F) -> F {
        F::from(a) + alpha * F::from(b) + alpha.square() * F::from(c)
    }

    fn s<F: PrimeField>(f: &Vec<F>, t: &Vec<F>) -> Vec<F> {
        let mut s = [f.clone(), t.clone()].concat();
        s.sort();
        s
    }

    // get difference vectors
    fn diff<F: PrimeField>(s: &Vec<F>) -> Vec<F> {
        (0..s.len() - 1).map(|i| s[i + 1] - s[i]).collect()
    }

    // check a ⊂ b
    fn multiset_check<F: PrimeField>(a: &Vec<F>, b: &Vec<F>) -> bool {
        a.iter().all(|vector: &F| b.contains(&vector))
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
        let alpha = Scalar::random(OsRng);
        let mut witness_vectors = witness_vectors(range, alpha);
        let xor_table = XORTable::<Scalar>::precompute();
        let mut t = xor_table.compress(alpha);

        // naive multset check
        assert!(multiset_check(&witness_vectors, &t));
        witness_vectors[0] += Scalar::one();
        assert!(!multiset_check(&witness_vectors, &t));

        witness_vectors[0] -= Scalar::one();
        t.sort();
        let s = s(&witness_vectors, &t);
        let s_prime = diff(&s);
        let t_prime = diff(&t);

        // f ⊂ t, s'
        assert!(multiset_check(&witness_vectors, &[t, s_prime].concat()));
    }
}
