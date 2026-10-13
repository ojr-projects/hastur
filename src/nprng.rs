//! numpy-compatible random numbers: `np.random.default_rng(seed)` = PCG64 (XSL-RR 128/64) seeded by SeedSequence,
//! with `random()`, `choice(k, size, p)` and `permutation(n)` reproduced draw for draw, so that the Rust port makes the
//! same synthetic Markov sets and the same permutations as proto/hlat.py.

const MULT: u128 = 0x2360ED051FC65DA44385DF649FCCF645;

pub struct Pcg64 {
    state: u128,
    inc: u128,
    has_u32: bool,
    u32buf: u32,
}

fn hashmix(value: u32, hash_const: &mut u32) -> u32 {
    let mut v = value ^ *hash_const;
    *hash_const = hash_const.wrapping_mul(0x931e8875);
    v = v.wrapping_mul(*hash_const);
    v ^ (v >> 16)
}

fn mix(x: u32, y: u32) -> u32 {
    let r = 0xca01f9ddu32.wrapping_mul(x).wrapping_sub(0x4973f715u32.wrapping_mul(y));
    r ^ (r >> 16)
}

/// numpy SeedSequence(entropy).generate_state(n_u64, uint64), for a non-negative integer entropy.
fn seed_sequence_u64(entropy: u64, n_u64: usize) -> Vec<u64> {
    let mut ent: Vec<u32> = Vec::new();
    let mut n = entropy;
    if n == 0 {
        ent.push(0);
    }
    while n > 0 {
        ent.push((n & 0xFFFF_FFFF) as u32);
        n >>= 32;
    }
    let mut pool = [0u32; 4];
    let mut hc = 0x43b0d7e5u32;
    for i in 0..4 {
        pool[i] = hashmix(if i < ent.len() { ent[i] } else { 0 }, &mut hc);
    }
    for i_src in 0..4 {
        for i_dst in 0..4 {
            if i_src != i_dst {
                let h = hashmix(pool[i_src], &mut hc);
                pool[i_dst] = mix(pool[i_dst], h);
            }
        }
    }
    for &e in ent.iter().skip(4) {
        for i_dst in 0..4 {
            let h = hashmix(e, &mut hc);
            pool[i_dst] = mix(pool[i_dst], h);
        }
    }
    let nw = 2 * n_u64;
    let mut st = vec![0u32; nw];
    let mut hb = 0x8b51f9ddu32;
    for i in 0..nw {
        let mut d = pool[i % 4] ^ hb;
        hb = hb.wrapping_mul(0x58f38ded);
        d = d.wrapping_mul(hb);
        d ^= d >> 16;
        st[i] = d;
    }
    (0..n_u64).map(|i| (st[2 * i] as u64) | ((st[2 * i + 1] as u64) << 32)).collect()
}

impl Pcg64 {
    pub fn new(seed: u64) -> Pcg64 {
        let v = seed_sequence_u64(seed, 4);
        let initstate = ((v[0] as u128) << 64) | v[1] as u128;
        let initseq = ((v[2] as u128) << 64) | v[3] as u128;
        let mut r = Pcg64 { state: 0, inc: (initseq << 1) | 1, has_u32: false, u32buf: 0 };
        r.step();
        r.state = r.state.wrapping_add(initstate);
        r.step();
        r
    }
    #[inline]
    fn step(&mut self) {
        self.state = self.state.wrapping_mul(MULT).wrapping_add(self.inc);
    }
    #[inline]
    pub fn next_u64(&mut self) -> u64 {
        self.step();
        let s = self.state;
        let x = ((s >> 64) as u64) ^ (s as u64);
        x.rotate_right((s >> 122) as u32)
    }
    #[inline]
    pub fn next_u32(&mut self) -> u32 {
        if self.has_u32 {
            self.has_u32 = false;
            return self.u32buf;
        }
        let n = self.next_u64();
        self.has_u32 = true;
        self.u32buf = (n >> 32) as u32;
        (n & 0xFFFF_FFFF) as u32
    }
    #[inline]
    pub fn random(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 * (1.0 / 9007199254740992.0)
    }
    /// Generator.choice(k, size=n, p=p): cdf = cumsum(p) / cdf[-1]; searchsorted(uniform, side='right').
    pub fn choice_p(&mut self, p: &[f64], n: usize) -> Vec<usize> {
        let mut cdf: Vec<f64> = Vec::with_capacity(p.len());
        let mut acc = 0.0;
        for &x in p {
            acc += x;
            cdf.push(acc);
        }
        let last = *cdf.last().unwrap();
        for c in cdf.iter_mut() {
            *c /= last;
        }
        (0..n).map(|_| {
            let u = self.random();
            cdf.iter().filter(|&&c| c <= u).count()
        }).collect()
    }
    fn interval(&mut self, max: u64) -> u64 {
        if max == 0 {
            return 0;
        }
        let mut mask = max;
        mask |= mask >> 1;
        mask |= mask >> 2;
        mask |= mask >> 4;
        mask |= mask >> 8;
        mask |= mask >> 16;
        mask |= mask >> 32;
        if max <= 0xFFFF_FFFF {
            loop {
                let v = (self.next_u32() as u64) & mask;
                if v <= max {
                    return v;
                }
            }
        } else {
            loop {
                let v = self.next_u64() & mask;
                if v <= max {
                    return v;
                }
            }
        }
    }
    /// Generator.permutation(n): arange(n) shuffled by numpy's Fisher-Yates (random_interval).
    pub fn permutation(&mut self, n: usize) -> Vec<usize> {
        let mut a: Vec<usize> = (0..n).collect();
        for i in (1..n).rev() {
            let j = self.interval(i as u64) as usize;
            a.swap(i, j);
        }
        a
    }
}

/// Reference values from numpy 2 (see rust/validate_full.py).
pub fn selftest() -> bool {
    let mut r = Pcg64::new(0);
    let a = [r.random(), r.random(), r.random()];
    let ok1 = a == [0.6369616873214543, 0.2697867137638703, 0.04097352393619469];
    let mut r = Pcg64::new(271828);
    let ok2 = r.permutation(10) == vec![4, 5, 1, 9, 3, 8, 6, 0, 2, 7];
    let mut r = Pcg64::new(0);
    let ok3 = r.choice_p(&[0.1, 0.2, 0.3, 0.4], 6) == vec![3, 1, 0, 0, 3, 3];
    eprintln!("nprng selftest: random {ok1}, permutation {ok2}, choice {ok3}");
    ok1 && ok2 && ok3
}
