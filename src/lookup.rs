//! plookup protocol
//! https://eprint.iacr.org/2020/315.pdf#page=6
use rand::rngs::OsRng;
use zkstd::common::{FftField, PrimeField};

mod table;

use table::XORTable;

pub(crate) struct Lookup<F: FftField> {
    a: Vec<F>,
    b: Vec<F>,
    c: Vec<F>,
}

impl<F: FftField> Lookup<F> {
    pub(crate) fn new(a: Vec<F>, b: Vec<F>, c: Vec<F>) -> Self {
        Self { a, b, c }
    }

    pub(crate) fn prove(&self, alpha: F, table: XORTable<F>) {
        let f = self.compress(alpha);
        let t = table.compress(alpha);
        let mut s = [f.clone(), t.clone()].concat();
        s.sort();
        let (h1, h2) = s.split_at(s.len() / 2);
        let (β, y) = (F::random(OsRng), F::random(OsRng));
        let z = compute_z(β, y, f, t, h1.to_vec(), h2.to_vec());
    }

    fn compress(&self, alpha: F) -> Vec<F> {
        let alpha2 = alpha.square();
        self.a
            .iter()
            .zip(self.b.iter())
            .zip(self.c.iter())
            .map(|((t1, t2), t3)| *t1 + alpha * t2 + alpha2 * t3)
            .collect()
    }
}

fn compute_z<F: FftField>(β: F, y: F, f: Vec<F>, t: Vec<F>, h1: Vec<F>, h2: Vec<F>) -> Vec<F> {
    let n = f.len();
    let one_β = F::one() + β;
    let mut z = vec![F::one()];
    let mut f_prev = F::one();
    let mut s_prev = F::one();
    for i in 2..n {
        let fi = compute_f(i, one_β, β, y, &f, &t);
        let gi = compute_g(i, one_β, β, y, &h1, &h2);
        f_prev *= fi;
        s_prev *= gi;
        let zi = f_prev / s_prev;
        z.push(zi)
    }
    z.push(F::one());
    z
}

fn compute_f<F: FftField>(i: usize, one_β: F, β: F, y: F, f: &Vec<F>, t: &Vec<F>) -> F {
    let left = y + f[i];
    let right = randomly_linear_combination(one_β, β, y, t[i], t[i + 1]);
    left * right
}

fn compute_g<F: FftField>(i: usize, one_β: F, β: F, y: F, h1: &Vec<F>, h2: &Vec<F>) -> F {
    let left = randomly_linear_combination(one_β, β, y, h1[i], h1[i + 1]);
    let right = randomly_linear_combination(one_β, β, y, h2[i], h2[i + 1]);
    left * right
}

// λ(1 + β) + a_i + β a_i_1
fn randomly_linear_combination<F: FftField>(one_β: F, β: F, y: F, a_i: F, a_i1: F) -> F {
    y * one_β + a_i + β * a_i1
}

#[cfg(test)]
mod tests {
    use super::table::XORTable;
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
        let t_2prime = [vec![Scalar::zero()], t_prime].concat();

        // s ⊂ t
        assert!(multiset_check(&s_prime, &t_2prime));
    }
}
