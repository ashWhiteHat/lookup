use rand::rngs::OsRng;
use zkstd::common::FftField;

pub(crate) fn challenge_scalar<F: FftField>() -> F {
    F::random(OsRng)
}
