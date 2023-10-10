use zkstd::common::PrimeField;

pub(crate) struct XORTable<F: PrimeField> {
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
        Self { a, b, c }
    }
}

#[cfg(test)]
mod tests {
    use super::XORTable;
    use bls_12_381::Fr as Scalar;

    #[test]
    fn table_generation() {
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
}
