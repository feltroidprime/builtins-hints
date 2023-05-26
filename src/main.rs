use lazy_static::lazy_static;
use num_bigint::{BigInt, RandBigInt, ToBigInt};
use num_integer::Integer;
use num_traits::{One, Signed, Zero};
use std::str::FromStr;

const BASE_BITS: u32 = 86;
const N_LIMBS: usize = 3;
const DEGREE: usize = N_LIMBS - 1;
const N_LIMBS_UNREDUCED: usize = 2 * N_LIMBS - 1;

lazy_static! {
    static ref BASE: BigInt = {
        // Perform the necessary computations to obtain the BigInt value
        BigInt::from(2u32).pow(BASE_BITS)
    };
    static ref P0:BigInt = {
        BigInt::from_str("69440356433466637143769089").unwrap()
    };

    static ref P1:BigInt = {
        BigInt::from_str("27625954992971143715037670").unwrap()
    };

    static ref P2:BigInt = {
        BigInt::from_str("3656382694611191768777988").unwrap()
    };

    static ref P:BigInt = {
        BigInt::from_str("21888242871839275222246405745257275088548364400416034343698204186575808495617").unwrap()
    };
    static ref P_LIMBS: Vec<BigInt> = {
        let mut limbs = Vec::with_capacity(N_LIMBS);
        limbs.push(P0.clone());
        limbs.push(P1.clone());
        limbs.push(P2.clone());
        limbs
    };

}

fn main() {
    println!("Hello, world!");
    println!("The value of BASE is: {}", *BASE);
    let base_2 = BASE.pow(2);
    println!("The value of BASE^2 is: {}", base_2);
    println!("The value of split(BASE²+4) is {:?}", split(base_2 + 4u32));
    let mut rng = rand::thread_rng();
    let a = rng.gen_bigint_range(&0.to_bigint().unwrap(), &P);
    let b = rng.gen_bigint_range(&0.to_bigint().unwrap(), &P);
    let a_limbs = split(a);
    let b_limbs = split(b);

    mul_hint(a_limbs, b_limbs);
    println!("ok")
}

// Inputs : 2 * N_LIMBS felts
// Ouputs : 2 * N_LIMBS felts + 2 * (UNREDUCED_N_LIMBS - 1) felts
fn mul_hint(
    a_limbs: Vec<BigInt>,
    b_limbs: Vec<BigInt>,
) -> (Vec<BigInt>, Vec<BigInt>, Vec<BigInt>, Vec<BigInt>) {
    let a = eval(&a_limbs);
    let b = eval(&b_limbs);
    let c = a * b;
    let (q, r) = c.div_rem(&P);
    let qs = split(q);
    let rs = split(r);
    let val_limbs = poly_mul(a_limbs, b_limbs);
    let q_p_plus_r_limbs = poly_mul_plus_c(&qs, &P_LIMBS.to_vec(), &rs);
    let diff_limbs = poly_sub(&q_p_plus_r_limbs, &val_limbs, N_LIMBS_UNREDUCED);
    let carries = reduce_zero_poly(&diff_limbs);
    let carries = abs_poly(&carries);
    let flag_poly = flag_poly(&diff_limbs);
    return (qs, rs, carries, flag_poly);
}

fn split(x: BigInt) -> Vec<BigInt> {
    let mut coeffs = Vec::with_capacity(N_LIMBS);
    let mut x = x;

    for n in (0..=DEGREE).rev() {
        let divisor = BASE.pow(n as u32);
        let (q, r) = x.div_rem(&divisor);
        coeffs.push(q);
        x = r;
    }
    coeffs.reverse();
    return coeffs;
}

fn eval(coeffs: &Vec<BigInt>) -> BigInt {
    let mut result = BigInt::zero();
    for i in 0..N_LIMBS {
        result += coeffs[i].clone() * BASE.pow(i as u32);
    }
    return result;
}

fn poly_mul(a: Vec<BigInt>, b: Vec<BigInt>) -> Vec<BigInt> {
    let mut result = vec![Zero::zero(); N_LIMBS_UNREDUCED];

    for i in 0..N_LIMBS {
        for j in 0..N_LIMBS {
            result[i + j] += a[i].clone() * b[j].clone();
        }
    }
    result
}

fn poly_mul_plus_c(a: &[BigInt], b: &[BigInt], c: &[BigInt]) -> Vec<BigInt> {
    let mut result = vec![BigInt::zero(); N_LIMBS_UNREDUCED];
    for i in 0..N_LIMBS {
        for j in 0..N_LIMBS {
            result[i + j] += a[i].clone() * b[j].clone();
        }
    }

    for i in 0..N_LIMBS {
        result[i] += c[i].clone();
    }
    return result;
}

fn poly_sub(a: &[BigInt], b: &[BigInt], n: usize) -> Vec<BigInt> {
    assert_eq!(a.len(), n);
    assert_eq!(b.len(), n);

    let mut result = vec![BigInt::zero(); n];
    for i in 0..n {
        result[i] = a[i].clone() - b[i].clone();
    }
    return result;
}

fn abs_poly(x: &[BigInt]) -> Vec<BigInt> {
    let mut result = Vec::with_capacity(x.len());

    for n in x {
        result.push(n.abs());
    }

    return result;
}

fn flag_poly(x: &[BigInt]) -> Vec<BigInt> {
    let mut result = Vec::with_capacity(N_LIMBS_UNREDUCED - 1);

    for n in &x[..N_LIMBS_UNREDUCED - 1] {
        if n >= &BigInt::zero() {
            result.push(BigInt::one());
        } else {
            result.push(BigInt::zero());
        }
    }

    result
}

fn reduce_zero_poly(x: &[BigInt]) -> Vec<BigInt> {
    let mut x = x.to_owned();
    let mut carries = vec![BigInt::zero(); x.len() - 1];

    for i in 0..x.len() - 1 {
        let (q, r) = x[i].clone().div_rem(&BASE);
        carries[i] = q;
        x[i] = r;
        assert!(x[i].is_zero());
        x[i + 1] += carries[i].clone();
    }
    assert!(x.last().unwrap().is_zero());
    return carries;
}
