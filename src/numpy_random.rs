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
//! The newer `Generator` (`numpy.random.default_rng`, the PCG64 bit
//! generator with the ziggurat normal distribution) is implemented by
//! [`NumpyGenerator`].

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

/// The multiplier of the PCG64 bit generator (numpy's
/// `PCG_DEFAULT_MULTIPLIER_128`).
const PCG64_MULTIPLIER: u128 = ((2549297995355413924u64 as u128) << 64) | 4865540595714422341u128;

/// The constants of numpy's `SeedSequence` (numpy `random/bit_generator.pyx`).
const INIT_A: u32 = 0x43b0_d7e5;
const MULT_A: u32 = 0x931e_8875;
const INIT_B: u32 = 0x8b51_f9dd;
const MULT_B: u32 = 0x58f3_8ded;
const MIX_MULT_L: u32 = 0xca01_f9dd;
const MIX_MULT_R: u32 = 0x4973_f715;
const XSHIFT: u32 = 16;

/// The size of the entropy pool of numpy's `SeedSequence`.
const POOL_SIZE: usize = 4;

/// numpy's `hashmix`.
fn hashmix(value: u32, hash_const: &mut u32) -> u32 {
    let mut value = value ^ *hash_const;
    *hash_const = hash_const.wrapping_mul(MULT_A);
    value = value.wrapping_mul(*hash_const);
    value ^ (value >> XSHIFT)
}

/// numpy's `mix`.
fn mix(x: u32, y: u32) -> u32 {
    let result = MIX_MULT_L
        .wrapping_mul(x)
        .wrapping_sub(MIX_MULT_R.wrapping_mul(y));
    result ^ (result >> XSHIFT)
}

/// The entropy pool of numpy's `SeedSequence` (the seed is broken up into
/// 32-bit words, the lowest bits first).
fn seed_sequence_pool(seed: u64) -> [u32; POOL_SIZE] {
    let mut entropy: Vec<u32> = Vec::new();
    let mut value = seed;
    while value > 0 {
        entropy.push((value & 0xffff_ffff) as u32);
        value >>= 32;
    }
    if entropy.is_empty() {
        entropy.push(0);
    }

    // Mix the entropy into the pool (numpy's `mix_entropy`)
    let mut mixer = [0u32; POOL_SIZE];
    let mut hash_const = INIT_A;
    for (i, value) in mixer.iter_mut().enumerate() {
        let entropy_value = entropy.get(i).copied().unwrap_or(0);
        *value = hashmix(entropy_value, &mut hash_const);
    }
    for i_src in 0..POOL_SIZE {
        let source = mixer[i_src];
        for (i_dst, value) in mixer.iter_mut().enumerate() {
            if i_src != i_dst {
                let hashed = hashmix(source, &mut hash_const);
                *value = mix(*value, hashed);
            }
        }
    }
    for &entropy_value in entropy.iter().skip(POOL_SIZE) {
        for value in mixer.iter_mut() {
            let hashed = hashmix(entropy_value, &mut hash_const);
            *value = mix(*value, hashed);
        }
    }

    mixer
}

/// The state of numpy's `SeedSequence.generate_state(n_words, uint64)`.
fn seed_sequence_state(seed: u64, n_words: usize) -> Vec<u64> {
    let pool = seed_sequence_pool(seed);
    let mut hash_const = INIT_B;
    let mut words: Vec<u32> = Vec::with_capacity(2 * n_words);
    for pool_value in pool.iter().cycle().take(2 * n_words) {
        let mut data_val = *pool_value;
        data_val ^= hash_const;
        hash_const = hash_const.wrapping_mul(MULT_B);
        data_val = data_val.wrapping_mul(hash_const);
        data_val ^= data_val >> XSHIFT;
        words.push(data_val);
    }

    // The 32-bit words are combined into 64-bit words (little endian)
    words
        .chunks(2)
        .map(|chunk| chunk[0] as u64 | ((chunk[1] as u64) << 32))
        .collect()
}

/// numpy's `Generator`: the PCG64 bit generator with the ziggurat normal
/// distribution (`numpy.random.default_rng`).
#[derive(Debug, Clone)]
pub struct NumpyGenerator {
    /// The 128-bit state of the bit generator
    state: u128,
    /// The 128-bit increment of the bit generator
    inc: u128,
}

impl NumpyGenerator {
    /// Create a generator with the given seed (numpy's
    /// `np.random.default_rng(seed)`).
    pub fn new(seed: u64) -> Self {
        // The state of the bit generator (numpy's `pcg64_set_seed`): the
        // first two words are the seed and the last two the increment; the
        // high halves come first
        let words = seed_sequence_state(seed, 4);
        let initstate = ((words[0] as u128) << 64) | words[1] as u128;
        let initseq = ((words[2] as u128) << 64) | words[3] as u128;

        // The increment (numpy's `pcg_setseq_128_srandom_r`)
        let inc = (((initseq >> 64) << 1 | (initseq & 0xffff_ffff_ffff_ffff) >> 63) << 64)
            | ((initseq & 0xffff_ffff_ffff_ffff) << 1 | 1);

        // The state
        let mut state: u128 = 0;
        state = state.wrapping_mul(PCG64_MULTIPLIER).wrapping_add(inc);
        state = state.wrapping_add(initstate);
        state = state.wrapping_mul(PCG64_MULTIPLIER).wrapping_add(inc);

        Self { state, inc }
    }

    /// The next 64-bit random number (numpy's `next_uint64`).
    pub fn next_u64(&mut self) -> u64 {
        self.state = self
            .state
            .wrapping_mul(PCG64_MULTIPLIER)
            .wrapping_add(self.inc);

        // The output function of the XSL-RR 128/64 generator
        let hi = (self.state >> 64) as u64;
        let lo = self.state as u64;
        let value = hi ^ lo;
        let rot = ((hi >> 58) & 63) as u32;
        value.rotate_right(rot)
    }

    /// A uniformly distributed number in [0, 1) (numpy's `random_double`).
    pub fn next_double(&mut self) -> Float {
        (self.next_u64() >> 11) as Float * (1.0 / 9_007_199_254_740_992.0)
    }

    /// A uniformly distributed number in [0, 1) (numpy's `random`).
    pub fn random_sample(&mut self) -> Float {
        self.next_double()
    }

    /// `n` uniformly distributed numbers in [0, 1).
    pub fn random_sample_vec(&mut self, n: usize) -> Vec<Float> {
        (0..n).map(|_| self.random_sample()).collect()
    }

    /// A normally distributed number (numpy's `standard_normal`: the
    /// ziggurat method).
    pub fn standard_normal(&mut self) -> Float {
        use crate::numpy_ziggurat::{FI, KI, NOR_INV_R, NOR_R, WI};

        loop {
            let mut r = self.next_u64();
            let idx = (r & 0xff) as usize;
            r >>= 8;
            let rabs = (r >> 1) & 0x000f_ffff_ffff_ffff;
            let x = rabs as Float * WI[idx] * if r & 1 != 0 { -1.0 } else { 1.0 };
            if rabs < KI[idx] {
                // 99.3 % of the time the value is returned here
                return x;
            }
            if idx == 0 {
                // The tail of the distribution
                loop {
                    let xx = -NOR_INV_R * (-self.next_double()).ln_1p();
                    let yy = -(-self.next_double()).ln_1p();
                    if yy + yy > xx * xx {
                        return if (rabs >> 8) & 1 != 0 {
                            -(NOR_R + xx)
                        } else {
                            NOR_R + xx
                        };
                    }
                }
            } else if (FI[idx - 1] - FI[idx]) * self.next_double() + FI[idx] < (-0.5 * x * x).exp()
            {
                return x;
            }
        }
    }

    /// `n` normally distributed numbers (numpy's `standard_normal`).
    pub fn standard_normal_vec(&mut self, n: usize) -> Vec<Float> {
        (0..n).map(|_| self.standard_normal()).collect()
    }

    /// A normally distributed number with the given mean and standard
    /// deviation (numpy's `normal`).
    pub fn normal(&mut self, loc: Float, scale: Float) -> Float {
        loc + scale * self.standard_normal()
    }

    /// A uniformly distributed number in [low, high) (numpy's `uniform`).
    pub fn uniform(&mut self, low: Float, high: Float) -> Float {
        low + (high - low) * self.random_sample()
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

    // The reference values of numpy's `Generator` are from:
    //   np.random.default_rng(0).standard_normal(8)
    //   np.random.default_rng(0).random(8)
    //   np.random.default_rng(42).standard_normal(4)
    //   np.random.default_rng(0).normal(0, 0.25, 6)
    #[test]
    fn test_generator_standard_normal_matches_numpy() {
        let mut rng = NumpyGenerator::new(0);
        let values = rng.standard_normal_vec(8);
        let expected = [
            0.1257302210933933,
            -0.1321048632913019,
            0.6404226504432821,
            0.10490011715303971,
            -0.535669373161111,
            0.36159505490948474,
            1.3040000451301372,
            0.9470809631292422,
        ];
        for (value, expected) in values.iter().zip(expected.iter()) {
            assert_relative_eq!(*value, *expected, epsilon = 1e-15);
        }
    }

    #[test]
    fn test_generator_random_matches_numpy() {
        let mut rng = NumpyGenerator::new(0);
        let values = rng.random_sample_vec(8);
        let expected = [
            0.6369616873214543,
            0.2697867137638703,
            0.04097352393619469,
            0.016527635528529094,
            0.8132702392002724,
            0.9127555772777217,
            0.6066357757671799,
            0.7294965609839984,
        ];
        for (value, expected) in values.iter().zip(expected.iter()) {
            assert_relative_eq!(*value, *expected, epsilon = 1e-15);
        }
    }

    #[test]
    fn test_generator_seeded_stream_matches_numpy() {
        let mut rng = NumpyGenerator::new(42);
        let values = rng.standard_normal_vec(4);
        let expected = [
            0.30471707975443135,
            -1.0399841062404955,
            0.7504511958064572,
            0.9405647163912139,
        ];
        for (value, expected) in values.iter().zip(expected.iter()) {
            assert_relative_eq!(*value, *expected, epsilon = 1e-15);
        }
    }

    #[test]
    fn test_generator_normal_matches_numpy() {
        let mut rng = NumpyGenerator::new(0);
        let values: Vec<Float> = (0..6).map(|_| rng.normal(0.0, 0.25)).collect();
        let expected = [
            0.031432555273348324,
            -0.03302621582282547,
            0.16010566261082051,
            0.026225029288259927,
            -0.13391734329027774,
            0.09039876372737118,
        ];
        for (value, expected) in values.iter().zip(expected.iter()) {
            assert_relative_eq!(*value, *expected, epsilon = 1e-15);
        }
    }
}
