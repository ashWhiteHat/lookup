use zkstd::common::FftField;

pub(crate) struct Polynomial<F: FftField> {
    coeffs: Vec<F>,
}

impl<F: FftField> Polynomial<F> {
    pub(crate) fn new(coeffs: Vec<F>) -> Self {
        Self { coeffs }
    }

    // log(n) inner product
    pub(crate) fn inner_product(&self, rhs: Self) -> F {
        assert_eq!(self.coeffs.len(), rhs.coeffs.len());
        self.coeffs
            .iter()
            .zip(rhs.coeffs.iter())
            .fold(F::zero(), |sum, (a, b)| sum + *a * *b)
    }
}
