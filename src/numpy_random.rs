//! A numpy-compatible random number generator.
//!
//! The FLORIS examples that use a random wind resource or a random layout
//! draw their numbers with numpy, so that the results are reproducible. This
//! module reimplements numpy's legacy `RandomState` (the MT19937 bit
//! generator with the polar method for the normal distribution), so that the
//! Rust examples draw exactly the same numbers as the Python examples.
//!
//! Reference: numpy `random/mtrand.pyx` and `random/distributions.c`
//! (`mt19937_seed`, `genrand_uint32`, `random_double`, `legacy_gauss`,
//! `buffered_bounded_masked_uint32`).
//!
//! Note: numpy's newer `Generator` (`numpy.random.default_rng`) uses the
//! PCG64 bit generator with the ziggurat normal distribution, which is a
//! different stream; it is not reproduced here.

use crate::types::Float;

const N: usize = 624;
const M: usize = 397;
const MATRIX_A: u32 = 0x9908_b0df;
const UPPER_MASK: u32 = 0x8000_0000;
const LOWER_MASK: u32 = 0x7fff_ffff;

/// The MT19937 bit generator (numpy's `mt19937`).
#[derive(Debug, Clone)]
pub struct Mt19937 {
    state: [u32; N],
    index: usize,
}

impl Mt19937 {
    /// Seed the generator (numpy's `mt19937_seed`, which calls `init_genrand`).
    pub fn new(seed: u32) -> Self {
        let mut generator = Self {
            state: [0; N],
            index: N,
        };
        generator.init_genrand(seed);
        generator
    }

    /// numpy's `init_genrand`.
    fn init_genrand(&mut self, seed: u32) {
        self.state[0] = seed;
        for i in 1..N {
            let previous = self.state[i - 1];
            self.state[i] = 1812433253u32
                .wrapping_mul(previous ^ (previous >> 30))
                .wrapping_add(i as u32);
        }
        self.index = N;
    }

    /// numpy's `genrand_uint32`.
    pub fn next_u32(&mut self) -> u32 {
        if self.index >= N {
            // The "twist": the new state is a linear combination of the
            // old state, shifted by M positions.
            let mut y: u32;
            for kk in 0..(N - M) {
                y = (self.state[kk] & UPPER_MASK) | (self.state[kk + 1] & LOWER_MASK);
                self.state[kk] =
                    self.state[kk + M] ^ (y >> 1) ^ if y & 1 != 0 { MATRIX_A } else { 0 };
            }
            for kk in (N - M)..(N - 1) {
                y = (self.state[kk] & UPPER_MASK) | (self.state[kk + 1] & LOWER_MASK);
                self.state[kk] =
                    self.state[kk + M - N] ^ (y >> 1) ^ if y & 1 != 0 { MATRIX_A } else { 0 };
            }
            y = (self.state[N - 1] & UPPER_MASK) | (self.state[0] & LOWER_MASK);
            self.state[N - 1] =
                self.state[M - 1] ^ (y >> 1) ^ if y & 1 != 0 { MATRIX_A } else { 0 };
            self.index = 0;
        }

        let mut y = self.state[self.index];
        self.index += 1;

        // numpy's tempering
        y ^= y >> 11;
        y ^= (y << 7) & 0x9d2c_5680;
        y ^= (y << 15) & 0xefc6_0000;
        y ^= y >> 18;
        y
    }

    /// numpy's `random_double`: two 32-bit draws give 53 random bits.
    pub fn next_double(&mut self) -> Float {
        let a = (self.next_u32() >> 5) as f64;
        let b = (self.next_u32() >> 6) as f64;
        (a * 67108864.0 + b) / 9007199254740992.0
    }
}

/// numpy's legacy `np.random.RandomState`.
#[derive(Debug, Clone)]
pub struct NumpyRandomState {
    bit_generator: Mt19937,
    /// The cached second value of the polar method (numpy's `has_gauss`).
    has_gauss: bool,
    gauss: Float,
}

impl NumpyRandomState {
    /// Create a random state with the given seed (numpy's `np.random.seed`).
    pub fn new(seed: u32) -> Self {
        Self {
            bit_generator: Mt19937::new(seed),
            has_gauss: false,
            gauss: 0.0,
        }
    }

    /// Reseed the generator (numpy's `np.random.seed`).
    pub fn seed(&mut self, seed: u32) {
        self.bit_generator = Mt19937::new(seed);
        self.has_gauss = false;
        self.gauss = 0.0;
    }

    /// A uniformly distributed number in [0, 1) (numpy's `random_sample`).
    pub fn random_sample(&mut self) -> Float {
        self.bit_generator.next_double()
    }

    /// A uniformly distributed number in [0, 1) (an alias of `random_sample`).
    pub fn random(&mut self) -> Float {
        self.random_sample()
    }

    /// `n` uniformly distributed numbers in [0, 1).
    pub fn random_sample_vec(&mut self, n: usize) -> Vec<Float> {
        (0..n).map(|_| self.random_sample()).collect()
    }

    /// A normally distributed number (numpy's `legacy_gauss`: the polar
    /// method, with the second value of the pair cached).
    pub fn standard_normal(&mut self) -> Float {
        if self.has_gauss {
            self.has_gauss = false;
            return self.gauss;
        }
        let (x1, x2, r2) = loop {
            let x1 = 2.0 * self.random_sample() - 1.0;
            let x2 = 2.0 * self.random_sample() - 1.0;
            let r2 = x1 * x1 + x2 * x2;
            if r2 < 1.0 && r2 != 0.0 {
                break (x1, x2, r2);
            }
        };
        let f = (-2.0 * r2.ln() / r2).sqrt();
        self.gauss = f * x1;
        self.has_gauss = true;
        f * x2
    }

    /// `n` normally distributed numbers (numpy's `np.random.randn`).
    pub fn standard_normal_vec(&mut self, n: usize) -> Vec<Float> {
        (0..n).map(|_| self.standard_normal()).collect()
    }

    /// A uniformly distributed number in [low, high) (numpy's `uniform`).
    pub fn uniform(&mut self, low: Float, high: Float) -> Float {
        low + (high - low) * self.random_sample()
    }

    /// `n` uniformly distributed numbers in [low, high).
    pub fn uniform_vec(&mut self, low: Float, high: Float, n: usize) -> Vec<Float> {
        (0..n).map(|_| self.uniform(low, high)).collect()
    }

    /// A random integer in [low, high) (numpy's `np.random.randint`).
    ///
    /// numpy draws 32-bit numbers and rejects the values that do not fit in
    /// the requested range (the "masked" method).
    pub fn randint(&mut self, low: i64, high: i64) -> i64 {
        assert!(high > low, "randint: high must be greater than low");
        let range = (high - low - 1) as u64;
        if range == 0 {
            return low;
        }
        // The smallest 2^k - 1 that is larger than or equal to the range
        let mut mask: u64 = 1;
        while mask < range {
            mask = (mask << 1) | 1;
        }
        let mask = mask as u32;
        let range = range as u32;
        loop {
            let value = self.bit_generator.next_u32();
            if value & mask <= range {
                return low + (value & mask) as i64;
            }
        }
    }

    /// A random element of `values`, drawn with the probabilities `p`
    /// (numpy's `np.random.choice` with `p`).
    pub fn choice(&mut self, values: &[Float], p: &[Float]) -> Float {
        assert_eq!(values.len(), p.len(), "choice: p must have the same length");
        let mut cdf: Vec<Float> = Vec::with_capacity(p.len());
        let mut total = 0.0;
        for probability in p {
            total += probability;
            cdf.push(total);
        }
        let last = *cdf.last().unwrap();
        for value in cdf.iter_mut() {
            *value /= last;
        }
        let sample = self.random_sample();
        // The first index whose cumulative probability is larger than the
        // sample (numpy's `searchsorted(..., side="right")`)
        let mut index = cdf.len() - 1;
        for (i, value) in cdf.iter().enumerate() {
            if *value > sample {
                index = i;
                break;
            }
        }
        values[index]
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    // The reference values are from numpy:
    //   np.random.seed(1); np.random.random_sample(5)
    //   np.random.seed(1); np.random.randn(5)
    //   np.random.seed(1); np.random.randint(0, 100, 5)
    //   np.random.seed(1); np.random.uniform(0, 1, 5)
    #[test]
    fn test_random_sample_matches_numpy() {
        let mut rng = NumpyRandomState::new(1);
        let values = rng.random_sample_vec(5);
        let expected = [
            0.417022004702574,
            0.7203244934421581,
            0.00011437481734488664,
            0.30233257263183977,
            0.14675589081711304,
        ];
        for (value, expected) in values.iter().zip(expected.iter()) {
            assert_relative_eq!(*value, *expected, epsilon = 1e-15);
        }
    }

    #[test]
    fn test_standard_normal_matches_numpy() {
        let mut rng = NumpyRandomState::new(1);
        let values = rng.standard_normal_vec(5);
        let expected = [
            1.6243453636632417,
            -0.6117564136500754,
            -0.5281717522634557,
            -1.0729686221561705,
            0.8654076293246785,
        ];
        for (value, expected) in values.iter().zip(expected.iter()) {
            assert_relative_eq!(*value, *expected, epsilon = 1e-15);
        }
    }

    #[test]
    fn test_randint_matches_numpy() {
        let mut rng = NumpyRandomState::new(1);
        let values: Vec<i64> = (0..5).map(|_| rng.randint(0, 100)).collect();
        assert_eq!(values, vec![37, 12, 72, 9, 75]);
    }

    #[test]
    fn test_uniform_matches_numpy() {
        let mut rng = NumpyRandomState::new(1);
        let values = rng.uniform_vec(0.0, 1.0, 5);
        let expected = [
            0.417022004702574,
            0.7203244934421581,
            0.00011437481734488664,
            0.30233257263183977,
            0.14675589081711304,
        ];
        for (value, expected) in values.iter().zip(expected.iter()) {
            assert_relative_eq!(*value, *expected, epsilon = 1e-15);
        }
    }

    #[test]
    fn test_choice_matches_numpy() {
        let mut rng = NumpyRandomState::new(1);
        let value = rng.choice(&[1.0, 2.0, 3.0], &[0.2, 0.3, 0.5]);
        assert_relative_eq!(value, 2.0, epsilon = 1e-15);
    }
}
