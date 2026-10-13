//! HASTUR, frozen F2 configuration (benchmarks/splits/freeze_v18.md). Ported from the Python prototype (not carried over),
//! with byte-identical outputs on every validation unit.
//!
//! What it computes, as F2 does:
//! - gapped k-of-w lattice (k = 4..5, w = 13) and contiguous L-mer table (L = 8), both strands, fg and bg;
//! - ranking statistics (v8): Wald double-centred z (gapped), soft-kernel Wald z (L-mer); genomic control on rank only;
//! - tested statistics: exact conditional binomial mid-p z (gapped), exact radius-1 Hamming-ball z (L-mer); no GC;
//! - seeds: top cells by the ranking statistic, de-duplicated as in Python;
//! - permutation max-T: R relabellings of the pooled sequences; only the permuted fg is counted (bg = pool - fg);
//!   Gumbel fit by moments; family thresholds at alpha 0.025.
//! - representation (with -o PREFIX; src/represent.rs): frame or EM PWM chosen on held-out data, then a site refit;
//!   writes PREFIX.tsv and PREFIX.meme.
//! Summary JSON on stdout.
//!
//! Counting: every strand of every set is flattened with N separators (>= w), so windows crossing a boundary hold an
//! N and drop out, as Python's padding does. Per position a 26-bit code of the 13-base window and a 13-bit N mask are
//! precomputed; a shape's cell is one PEXT of the window code, so a shape costs one pass of mask-test + PEXT +
//! increment per window.
//!
//! Usage: hastur <fg.fa> <bg.fa> [-o PREFIX] [--perm R] [--seed S] [--threads T] [--dump-cells N]

mod nprng;
mod represent;

use rayon::prelude::*;
use std::io::{BufRead, BufReader};
use std::time::Instant;

const W: usize = 13;
const KMIN: usize = 4;
const KMAX: usize = 5;
const L: usize = 8;
const SOFT_B: f64 = 0.05;
const N_SEEDS: usize = 12;
const SEP: usize = 16;
const ALPHA: f64 = 0.025;

// ---------------------------------------------------------------- input
fn read_fasta(path: &str) -> Vec<Vec<u8>> {
    let f = std::fs::File::open(path).unwrap_or_else(|e| panic!("{path}: {e}"));
    let mut out: Vec<Vec<u8>> = Vec::new();
    let mut cur: Option<Vec<u8>> = None;
    for line in BufReader::new(f).lines() {
        let line = line.unwrap();
        if line.starts_with('>') {
            if let Some(s) = cur.take() {
                out.push(s);
            }
            cur = Some(Vec::new());
        } else if let Some(s) = cur.as_mut() {
            s.extend(line.trim().bytes().map(|b| match b {
                b'A' | b'a' => 0u8,
                b'C' | b'c' => 1,
                b'G' | b'g' => 2,
                b'T' | b't' => 3,
                _ => 4,
            }));
        }
    }
    if let Some(s) = cur.take() {
        out.push(s);
    }
    out
}

fn revcomp(s: &[u8]) -> Vec<u8> {
    s.iter().rev().map(|&c| if c < 4 { 3 - c } else { 4 }).collect()
}

/// Window codes of a flattened strand: (26-bit code of bases t..t+12, 13-bit N mask; bit 12-j = base t+j is N).
struct Strand {
    code: Vec<u32>,
    nmask: Vec<u16>,
}

fn strand(seqs: &[&[u8]]) -> Strand {
    let mut a: Vec<u8> = Vec::with_capacity(seqs.iter().map(|s| s.len() + SEP).sum::<usize>() + W);
    for s in seqs {
        a.extend_from_slice(s);
        a.extend(std::iter::repeat(4u8).take(SEP));
    }
    let n = a.len();
    let mut code = vec![0u32; n];
    let mut nmask = vec![0u16; n];
    let (mut c, mut m) = (0u32, 0u16);
    let get = |i: usize| if i < n { a[i] } else { 4 };
    for j in 0..W {
        let b = get(j);
        c = (c << 2) | (b & 3) as u32;
        m = (m << 1) | (b == 4) as u16;
    }
    let cmask: u32 = (1u32 << (2 * W)) - 1;
    let mmask: u16 = (1u16 << W) - 1;
    for t in 0..n {
        code[t] = c;
        nmask[t] = m;
        let b = get(t + W);
        c = ((c << 2) | (b & 3) as u32) & cmask;
        m = ((m << 1) | (b == 4) as u16) & mmask;
    }
    Strand { code, nmask }
}

/// Both strands of a set of sequences.
struct SetS {
    fwd: Strand,
    rev: Strand,
}

fn set_of(seqs: &[&[u8]]) -> SetS {
    let rc: Vec<Vec<u8>> = seqs.iter().map(|s| revcomp(s)).collect();
    let rcr: Vec<&[u8]> = rc.iter().map(|v| v.as_slice()).collect();
    SetS { fwd: strand(seqs), rev: strand(&rcr) }
}

// ---------------------------------------------------------------- counting
#[inline]
fn pext(x: u32, m: u32) -> u32 {
    #[cfg(all(target_arch = "x86_64", target_feature = "bmi2"))]
    unsafe {
        return core::arch::x86_64::_pext_u32(x, m);
    }
    #[allow(unreachable_code)]
    {
        let (mut r, mut bb, mut mm) = (0u32, 1u32, m);
        while mm != 0 {
            let low = mm & mm.wrapping_neg();
            if x & low != 0 {
                r |= bb;
            }
            bb <<= 1;
            mm &= mm - 1;
        }
        r
    }
}

#[derive(Clone)]
struct Shape {
    offs: Vec<usize>,
    k: usize,
    bm: u32, // 2-bit selection mask over the 26-bit window code (offset o -> bits 2*(12-o), +1)
    pm: u16, // N-mask bits of the shape's positions
}

fn shape(offs: &[usize]) -> Shape {
    let mut bm = 0u32;
    let mut pm = 0u16;
    for &o in offs {
        bm |= 3u32 << (2 * (W - 1 - o));
        pm |= 1u16 << (W - 1 - o);
    }
    Shape { offs: offs.to_vec(), k: offs.len(), bm, pm }
}

fn all_shapes(k: usize) -> Vec<Shape> {
    // (0,) + combinations(range(1, W), k - 1), in itertools order
    fn rec(start: usize, left: usize, cur: &mut Vec<usize>, out: &mut Vec<Vec<usize>>) {
        if left == 0 {
            out.push(cur.clone());
            return;
        }
        for v in start..W {
            cur.push(v);
            rec(v + 1, left - 1, cur, out);
            cur.pop();
        }
    }
    let mut out = Vec::new();
    rec(1, k - 1, &mut vec![0], &mut out);
    out.iter().map(|o| shape(o)).collect()
}

fn count_one(st: &Strand, sh: &Shape, hist: &mut [u32]) {
    let (bm, pm) = (sh.bm, sh.pm);
    for (c, m) in st.code.iter().zip(st.nmask.iter()) {
        if m & pm == 0 {
            hist[pext(*c, bm) as usize] += 1;
        }
    }
}

/// Counts of every cell of every shape, both strands: Vec over shapes of Vec<f64> (4^k).
fn count_shapes(s: &SetS, shapes: &[Shape]) -> Vec<Vec<f64>> {
    shapes
        .par_iter()
        .map(|sh| {
            let mut h = vec![0u32; 1 << (2 * sh.k)];
            count_one(&s.fwd, sh, &mut h);
            count_one(&s.rev, sh, &mut h);
            h.into_iter().map(|v| v as f64).collect()
        })
        .collect()
}

// ---------------------------------------------------------------- Markov order 1
#[derive(Clone)]
struct Model {
    pi: [f64; 4],
    t: [[f64; 4]; 4],
}

fn markov_from_di(di: &[f64]) -> Model {
    let mut c = [[0.0; 4]; 4];
    for a in 0..4 {
        for b in 0..4 {
            c[a][b] = di[a * 4 + b] + 0.5;
        }
    }
    let tot: f64 = c.iter().flatten().sum();
    let mut pi = [0.0; 4];
    let mut t = [[0.0; 4]; 4];
    for a in 0..4 {
        let rs: f64 = c[a].iter().sum();
        pi[a] = rs / tot;
        for b in 0..4 {
            t[a][b] = c[a][b] / rs;
        }
    }
    Model { pi, t }
}

fn matmul(a: &[[f64; 4]; 4], b: &[[f64; 4]; 4]) -> [[f64; 4]; 4] {
    let mut r = [[0.0; 4]; 4];
    for i in 0..4 {
        for j in 0..4 {
            r[i][j] = (0..4).map(|k| a[i][k] * b[k][j]).sum();
        }
    }
    r
}

/// P(pattern) for every letter assignment of a shape under the order-1 chain (closed form, Python's index order).
fn expect_probs(offs: &[usize], m: &Model, tpow: &[[[f64; 4]; 4]]) -> Vec<f64> {
    let mut v: Vec<f64> = m.pi.to_vec();
    for i in 1..offs.len() {
        let g = &tpow[offs[i] - offs[i - 1]];
        let mut nv = vec![0.0; v.len() * 4];
        for (idx, &x) in v.iter().enumerate() {
            let last = idx & 3;
            for b in 0..4 {
                nv[idx * 4 + b] = x * g[last][b];
            }
        }
        v = nv;
    }
    v
}

fn tpowers(m: &Model, upto: usize) -> Vec<[[f64; 4]; 4]> {
    let mut p = vec![[[0.0; 4]; 4]; upto + 1];
    for i in 0..4 {
        p[0][i][i] = 1.0;
    }
    for g in 1..=upto {
        p[g] = matmul(&p[g - 1], &m.t);
    }
    p
}

// ---------------------------------------------------------------- special functions
fn ln_gamma(x: f64) -> f64 {
    // Lanczos (g = 7, n = 9)
    const C: [f64; 9] = [
        0.99999999999980993, 676.5203681218851, -1259.1392167224028, 771.32342877765313, -176.61502916214059,
        12.507343278686905, -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7,
    ];
    if x < 0.5 {
        return (std::f64::consts::PI / (std::f64::consts::PI * x).sin()).ln() - ln_gamma(1.0 - x);
    }
    let x = x - 1.0;
    let mut a = C[0];
    let t = x + 7.5;
    for (i, c) in C.iter().enumerate().skip(1) {
        a += c / (x + i as f64);
    }
    0.5 * (2.0 * std::f64::consts::PI).ln() + (x + 0.5) * t.ln() - t + a.ln()
}

fn betacf(a: f64, b: f64, x: f64) -> f64 {
    let (fpmin, eps) = (1e-300, 1e-15);
    let (qab, qap, qam) = (a + b, a + 1.0, a - 1.0);
    let mut c = 1.0;
    let mut d = 1.0 - qab * x / qap;
    if d.abs() < fpmin {
        d = fpmin;
    }
    d = 1.0 / d;
    let mut h = d;
    for m in 1..10000 {
        let m = m as f64;
        let m2 = 2.0 * m;
        let aa = m * (b - m) * x / ((qam + m2) * (a + m2));
        d = 1.0 + aa * d;
        if d.abs() < fpmin {
            d = fpmin;
        }
        c = 1.0 + aa / c;
        if c.abs() < fpmin {
            c = fpmin;
        }
        d = 1.0 / d;
        h *= d * c;
        let aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2));
        d = 1.0 + aa * d;
        if d.abs() < fpmin {
            d = fpmin;
        }
        c = 1.0 + aa / c;
        if c.abs() < fpmin {
            c = fpmin;
        }
        d = 1.0 / d;
        let del = d * c;
        h *= del;
        if (del - 1.0).abs() < eps {
            break;
        }
    }
    h
}

/// Regularised incomplete beta I_x(a, b).
fn betai(a: f64, b: f64, x: f64) -> f64 {
    if x <= 0.0 {
        return 0.0;
    }
    if x >= 1.0 {
        return 1.0;
    }
    let lbt = ln_gamma(a + b) - ln_gamma(a) - ln_gamma(b) + a * x.ln() + b * (1.0 - x).ln();
    if x < (a + 1.0) / (a + b + 2.0) {
        (lbt.exp() * betacf(a, b, x) / a).min(1.0)
    } else {
        (1.0 - lbt.exp() * betacf(b, a, 1.0 - x) / b).max(0.0)
    }
}

/// Mid-p upper tail of Binomial(n, p) at f: P(X > f) + P(X = f) / 2 (scipy: binom.sf(f) + 0.5 binom.pmf(f)).
fn binom_midp(f: f64, n: f64, p: f64) -> f64 {
    if n <= 0.0 {
        return 0.5; // scipy: sf(0; 0, p) + 0.5 pmf(0; 0, p) = 0 + 0.5
    }
    let gt = if f >= n { 0.0 } else { betai(f + 1.0, n - f, p) };
    let lpmf = ln_gamma(n + 1.0) - ln_gamma(f + 1.0) - ln_gamma(n - f + 1.0)
        + if f > 0.0 { f * p.ln() } else { 0.0 }
        + if n - f > 0.0 { (n - f) * (1.0 - p).ln() } else { 0.0 };
    gt + 0.5 * lpmf.exp()
}

/// Upper-tail inverse normal z = Phi^-1(1 - p) (Acklam), p clipped to [1e-300, 1] as Python does.
fn isf(p: f64) -> f64 {
    let p = p.clamp(1e-300, 1.0);
    let q = 1.0 - p; // lower-tail probability of -z ... use symmetric form on the smaller tail
    let x = if p < 0.5 { p } else { q };
    const A: [f64; 6] = [-3.969683028665376e+01, 2.209460984245205e+02, -2.759285104469687e+02, 1.383577518672690e+02,
        -3.066479806614716e+01, 2.506628277459239e+00];
    const B: [f64; 5] = [-5.447609879822406e+01, 1.615858368580409e+02, -1.556989798598866e+02, 6.680131188771972e+01,
        -1.328068155288572e+01];
    const C: [f64; 6] = [-7.784894002430293e-03, -3.223964580411365e-01, -2.400758277161838e+00, -2.549732539343734e+00,
        4.374664141464968e+00, 2.938163982698783e+00];
    const D: [f64; 4] = [7.784695709041462e-03, 3.224671290700398e-01, 2.445134137142996e+00, 3.754408661907416e+00];
    // lower-tail quantile of x (x <= 0.5) -> negative value
    let lq = if x < 0.02425 {
        let r = (-2.0 * x.ln()).sqrt();
        (((((C[0] * r + C[1]) * r + C[2]) * r + C[3]) * r + C[4]) * r + C[5]) / ((((D[0] * r + D[1]) * r + D[2]) * r + D[3]) * r + 1.0)
    } else {
        let r = x - 0.5;
        let s = r * r;
        (((((A[0] * s + A[1]) * s + A[2]) * s + A[3]) * s + A[4]) * s + A[5]) * r
            / (((((B[0] * s + B[1]) * s + B[2]) * s + B[3]) * s + B[4]) * s + 1.0)
    };
    // z = Phi^-1(1 - p): if p < 0.5, z = -lq(p) > 0; else z = lq(1 - p) (<= 0)
    if p < 0.5 { -lq } else { lq }
}

/// erfc (Numerical Recipes erfcc, rel. error < 1.2e-7): only for the Bonferroni fallback of p_adj.
fn libm_erfc(x: f64) -> f64 {
    let z = x.abs();
    let t = 1.0 / (1.0 + 0.5 * z);
    let ans = t * (-z * z - 1.26551223 + t * (1.00002368 + t * (0.37409196 + t * (0.09678418 + t * (-0.18628806
        + t * (0.27886807 + t * (-1.13520398 + t * (1.48851587 + t * (-0.82215223 + t * 0.17087277))))))))).exp();
    if x >= 0.0 { ans } else { 2.0 - ans }
}

// ---------------------------------------------------------------- statistics
#[inline]
fn exact_z(f: f64, n: f64, p0: f64, pal: bool) -> f64 {
    let z = isf(binom_midp(f.round(), n.round(), p0));
    if pal { z / std::f64::consts::SQRT_2 } else { z }
}

/// Upper bound on exact_z: mid-p >= pmf(f) / 2, so z <= isf(pmf / 2). Three log-gammas instead of a continued fraction.
#[inline]
fn exact_z_ub(f: f64, n: f64, p0: f64, pal: bool) -> f64 {
    let (f, n) = (f.round(), n.round());
    if n <= 0.0 {
        return 0.0;
    }
    let lpmf = ln_gamma(n + 1.0) - ln_gamma(f + 1.0) - ln_gamma(n - f + 1.0)
        + if f > 0.0 { f * p0.ln() } else { 0.0 }
        + if n - f > 0.0 { (n - f) * (1.0 - p0).ln() } else { 0.0 };
    let z = isf(0.5 * lpmf.exp());
    if pal { z / std::f64::consts::SQRT_2 } else { z }
}

/// Exact maximum of exact_z over cells (f, n, p0, pal), evaluating the exact tail only while a cell's upper bound can
/// still beat the best value found. Identical to the full maximum (rigorous bound), far fewer tail evaluations.
fn max_exact_z(cells: &[(f64, f64, f64, bool)]) -> f64 {
    let mut ub: Vec<(f64, usize)> = cells.iter().enumerate().map(|(i, c)| (exact_z_ub(c.0, c.1, c.2, c.3), i)).collect();
    ub.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap_or(std::cmp::Ordering::Equal));
    let mut best = f64::NEG_INFINITY;
    for (u, i) in ub {
        if u <= best {
            break;
        }
        let c = cells[i];
        let z = exact_z(c.0, c.1, c.2, c.3);
        if z > best {
            best = z;
        }
    }
    best
}

fn palindromic(offs: &[usize]) -> Option<Vec<bool>> {
    let span = *offs.last().unwrap();
    let mut mir: Vec<usize> = offs.iter().map(|o| span - o).collect();
    mir.sort();
    if mir != offs {
        return None;
    }
    let k = offs.len();
    Some((0..(1usize << (2 * k))).map(|i| {
        let d: Vec<usize> = (0..k).map(|j| (i >> (2 * (k - 1 - j))) & 3).collect();
        (0..k).all(|j| d[j] == 3 - d[k - 1 - j])
    }).collect())
}

struct GappedOut {
    zrank: Vec<f64>,
    ok: Vec<bool>,
    max_test: f64,
    cells: Vec<(f64, f64, f64, bool)>,
}

/// Rank (Wald) double-centred z for one shape, and the exact maximum of the tested (exact mid-p) z.
fn gapped_stats(f: &[f64], b: &[f64], pf: &[f64], pb: &[f64], pal: &Option<Vec<bool>>, want_rank: bool) -> GappedOut {
    let nf: f64 = f.iter().sum();
    let nb: f64 = b.iter().sum();
    let n = f.len();
    let mut zr = vec![0.0; if want_rank { n } else { 0 }];
    let mut ok = vec![false; if want_rank { n } else { 0 }];
    let mut cells = Vec::with_capacity(n);
    let s2 = std::f64::consts::SQRT_2;
    for i in 0..n {
        let ef = nf * pf[i] + 1e-9;
        let eb = nb * pb[i] + 1e-9;
        let pal_i = pal.as_ref().map_or(false, |m| m[i]);
        cells.push((f[i], f[i] + b[i], ef / (ef + eb), pal_i));
        if want_rank {
            let d = ((f[i] + 0.5) / ef).ln() - ((b[i] + 0.5) / eb).ln();
            let mut z = d / (1.0 / (f[i] + 0.5) + 1.0 / (b[i] + 0.5)).sqrt();
            if pal_i {
                z /= s2;
            }
            zr[i] = z;
            ok[i] = nf * pf[i] > 5.0 && nb * pb[i] > 5.0;
        }
    }
    let max_test = max_exact_z(&cells);
    GappedOut { zrank: zr, ok, max_test, cells }
}

fn mode_smooth(t: &[f64], k: &[[f64; 4]; 4]) -> Vec<f64> {
    let mut cur = t.to_vec();
    for j in 0..L {
        let stride = 1usize << (2 * (L - 1 - j));
        let mut nxt = vec![0.0; cur.len()];
        for idx in 0..cur.len() {
            let a = (idx / stride) & 3;
            let base = idx - a * stride;
            let mut s = 0.0;
            for b in 0..4 {
                s += k[a][b] * cur[base + b * stride];
            }
            nxt[idx] = s;
        }
        cur = nxt;
    }
    cur
}

fn ball(v: &[f64]) -> Vec<f64> {
    let mut out: Vec<f64> = v.iter().map(|x| -((L - 1) as f64) * x).collect();
    for j in 0..L {
        let stride = 1usize << (2 * (L - 1 - j));
        for idx in 0..v.len() {
            let a = (idx / stride) & 3;
            let base = idx - a * stride;
            out[idx] += (0..4).map(|b| v[base + b * stride]).sum::<f64>();
        }
    }
    out
}

struct LmerOut {
    rank: Vec<f64>,
    cells: Vec<(f64, f64, f64, bool)>,
    ok: Vec<bool>,
    max_test: f64,
}

fn lmer_stats(fl: &[f64], bl: &[f64], mf: &Model, mb: &Model, want_rank: bool) -> LmerOut {
    let offs: Vec<usize> = (0..L).collect();
    let pf = expect_probs(&offs, mf, &tpowers(mf, 1));
    let pb = expect_probs(&offs, mb, &tpowers(mb, 1));
    let (nf, nb): (f64, f64) = (fl.iter().sum(), bl.iter().sum());
    let muf: Vec<f64> = pf.iter().map(|p| p * nf).collect();
    let mub: Vec<f64> = pb.iter().map(|p| p * nb).collect();
    let pal = palindromic(&offs).unwrap();
    let s2 = std::f64::consts::SQRT_2;
    let (bf, bb, bmf, bmb) = (ball(fl), ball(bl), ball(&muf), ball(&mub));
    let cells: Vec<(f64, f64, f64, bool)> = (0..fl.len()).map(|i| {
        (bf[i], bf[i] + bb[i], (bmf[i] + 1e-12) / (bmf[i] + bmb[i] + 2e-12), pal[i])
    }).collect();
    let max_test = {
        // parallel pruned maximum over chunks
        cells.par_chunks(4096).map(max_exact_z).reduce(|| f64::NEG_INFINITY, f64::max)
    };
    if !want_rank {
        return LmerOut { rank: vec![], cells, ok: vec![], max_test };
    }
    let k1 = {
        let mut k = [[SOFT_B; 4]; 4];
        for i in 0..4 {
            k[i][i] = 1.0 - SOFT_B + SOFT_B;
        }
        k
    };
    let mut k2 = k1;
    for i in 0..4 {
        for j in 0..4 {
            k2[i][j] = k1[i][j] * k1[i][j];
        }
    }
    let (sf, sb) = (mode_smooth(fl, &k1), mode_smooth(bl, &k1));
    let (ef, eb) = (mode_smooth(&muf, &k1), mode_smooth(&mub, &k1));
    let (vf, vb) = (mode_smooth(fl, &k2), mode_smooth(bl, &k2));
    let mut rank = vec![0.0; fl.len()];
    let mut ok = vec![false; fl.len()];
    for i in 0..fl.len() {
        let d = ((sf[i] + 0.5) / (ef[i] + 1e-9)).ln() - ((sb[i] + 0.5) / (eb[i] + 1e-9)).ln();
        let mut z = d / (vf[i] / (sf[i] + 0.5).powi(2) + vb[i] / (sb[i] + 0.5).powi(2) + 1e-12).sqrt();
        if pal[i] {
            z /= s2;
        }
        rank[i] = z;
        ok[i] = ef[i] > 5.0 && eb[i] > 5.0;
    }
    LmerOut { rank, cells, ok, max_test }
}

fn median(v: &mut Vec<f64>) -> f64 {
    if v.is_empty() {
        return 0.0;
    }
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let n = v.len();
    if n % 2 == 1 { v[n / 2] } else { 0.5 * (v[n / 2 - 1] + v[n / 2]) }
}

// ---------------------------------------------------------------- patterns (Python's pattern_str / consistent_overlap)
fn pattern_str(offs: &[usize], cell: usize) -> String {
    let k = offs.len();
    let mut s = vec![b'N'; offs[k - 1] + 1];
    for (j, &o) in offs.iter().enumerate() {
        s[o] = b"ACGT"[(cell >> (2 * (k - 1 - j))) & 3];
    }
    String::from_utf8(s).unwrap()
}

fn consistent_overlap(p: &str, q: &str, min_match: usize) -> bool {
    let p = p.as_bytes();
    let rc: Vec<u8> = q.bytes().rev().map(|c| match c { b'A' => b'T', b'C' => b'G', b'G' => b'C', b'T' => b'A', x => x }).collect();
    let cnt = |s: &[u8]| s.iter().filter(|&&c| c != b'N').count();
    for c in [q.as_bytes().to_vec(), rc] {
        let need = min_match.min(cnt(&c)).min(cnt(p));
        for off in -(c.len() as isize) + 1..p.len() as isize {
            let (mut agree, mut conflict) = (0usize, 0usize);
            for (j, &ch) in c.iter().enumerate() {
                let i = off + j as isize;
                if i >= 0 && (i as usize) < p.len() && ch != b'N' && p[i as usize] != b'N' {
                    if ch == p[i as usize] { agree += 1 } else { conflict += 1 }
                }
            }
            if conflict == 0 && agree >= need {
                return true;
            }
        }
    }
    false
}

// ---------------------------------------------------------------- the detection pass
struct Pass {
    gapped_max_test: f64,
    lmer_max_test: f64,
    seeds: Vec<(String, f64, f64, String)>, // (family, z_rank, z_test, pattern)
    lam_rank_g: f64,
    lam_rank_l: f64,
}

struct Counts {
    g: Vec<Vec<Vec<f64>>>, // per k (index k - KMIN), per shape, cells
    l: Vec<f64>,
    di: Vec<f64>,
}

fn count_all(s: &SetS, shapes: &[Vec<Shape>], lsh: &Shape, dsh: &Shape) -> Counts {
    let g = shapes.iter().map(|sh| count_shapes(s, sh)).collect();
    let l = count_shapes(s, std::slice::from_ref(lsh)).pop().unwrap();
    let di = count_shapes(s, std::slice::from_ref(dsh)).pop().unwrap();
    Counts { g, l, di }
}

fn sub(a: &Counts, b: &Counts) -> Counts {
    let s = |x: &Vec<f64>, y: &Vec<f64>| x.iter().zip(y).map(|(u, v)| u - v).collect::<Vec<f64>>();
    Counts {
        g: a.g.iter().zip(&b.g).map(|(ka, kb)| ka.iter().zip(kb).map(|(u, v)| s(u, v)).collect()).collect(),
        l: s(&a.l, &b.l),
        di: s(&a.di, &b.di),
    }
}

fn detect(f: &Counts, b: &Counts, shapes: &[Vec<Shape>], want_seeds: bool) -> Pass {
    let (mf, mb) = (markov_from_di(&f.di), markov_from_di(&b.di));
    let (tf, tb) = (tpowers(&mf, W), tpowers(&mb, W));
    // gapped: per shape, in parallel
    let per: Vec<(f64, Vec<(f64, usize, usize, usize, f64)>, Vec<f64>)> = shapes
        .iter()
        .enumerate()
        .flat_map(|(ki, sh)| sh.iter().enumerate().map(move |(si, s)| (ki, si, s)))
        .collect::<Vec<_>>()
        .par_iter()
        .map(|&(ki, si, s)| {
            let pf = expect_probs(&s.offs, &mf, &tf);
            let pb = expect_probs(&s.offs, &mb, &tb);
            let pal = palindromic(&s.offs);
            let o = gapped_stats(&f.g[ki][si], &b.g[ki][si], &pf, &pb, &pal, want_seeds);
            let mx = o.max_test;
            if !want_seeds {
                return (mx, vec![], vec![]);
            }
            let mut idx: Vec<usize> = (0..o.zrank.len()).collect();
            idx.sort_by(|&a, &b| o.zrank[b].partial_cmp(&o.zrank[a]).unwrap());
            let cand = idx.iter().take(40).map(|&c| {
                let x = o.cells[c];
                (o.zrank[c], ki, si, c, exact_z(x.0, x.1, x.2, x.3))
            }).collect();
            let zok: Vec<f64> = o.zrank.iter().zip(&o.ok).filter(|(_, &k)| k).map(|(z, _)| *z).collect();
            (mx, cand, zok)
        })
        .collect();
    let gmax = per.iter().map(|p| p.0).fold(f64::NEG_INFINITY, f64::max);
    let lo = lmer_stats(&f.l, &b.l, &mf, &mb, want_seeds);
    let lmax = lo.max_test;
    if !want_seeds {
        return Pass { gapped_max_test: gmax, lmer_max_test: lmax, seeds: vec![], lam_rank_g: 1.0, lam_rank_l: 1.0 };
    }
    let mut zall: Vec<f64> = per.iter().flat_map(|p| p.2.iter().map(|z| z * z)).collect();
    let lam_g = (median(&mut zall) / 0.4549).max(1.0);
    let mut cand: Vec<(f64, usize, usize, usize, f64)> = per.into_iter().flat_map(|p| p.1).collect();
    cand.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap());
    let mut seeds: Vec<(String, f64, f64, String)> = Vec::new();
    for (z, ki, si, c, zt) in cand {
        let s = pattern_str(&shapes[ki][si].offs, c);
        if seeds.iter().any(|t| consistent_overlap(&s, &t.3, 4)) {
            continue;
        }
        seeds.push(("gapped".into(), z / lam_g.sqrt(), zt, s));
        if seeds.len() >= N_SEEDS {
            break;
        }
    }
    let mut zl: Vec<f64> = lo.rank.iter().zip(&lo.ok).filter(|(_, &k)| k).map(|(z, _)| z * z).collect();
    let lam_l = if zl.is_empty() { 1.0 } else { (median(&mut zl) / 0.4549).max(1.0) };
    let mut idx: Vec<usize> = (0..lo.rank.len()).collect();
    idx.sort_by(|&a, &b| lo.rank[b].partial_cmp(&lo.rank[a]).unwrap());
    let loffs: Vec<usize> = (0..L).collect();
    let mut nl = 0;
    for &i in idx.iter().take(2000) {
        let s = pattern_str(&loffs, i);
        if seeds.iter().any(|t| t.0 == "lmer" && consistent_overlap(&s, &t.3, 5)) {
            continue;
        }
        let x = lo.cells[i];
        seeds.push(("lmer".into(), lo.rank[i] / lam_l.sqrt(), exact_z(x.0, x.1, x.2, x.3), s));
        nl += 1;
        if nl >= N_SEEDS {
            break;
        }
    }
    Pass { gapped_max_test: gmax, lmer_max_test: lmax, seeds, lam_rank_g: lam_g, lam_rank_l: lam_l }
}

/// Gumbel fit by moments (numpy mean / std(ddof=1) with pairwise sums): (mu, beta, threshold at ALPHA).
fn gumbel(m: &[f64]) -> (f64, f64, f64) {
    let n = m.len() as f64;
    let mean = represent::np_sum(m) / n;
    let sq: Vec<f64> = m.iter().map(|x| (x - mean) * (x - mean)).collect();
    let sd = (represent::np_sum(&sq) / (n - 1.0)).sqrt();
    let beta = sd * 6f64.sqrt() / std::f64::consts::PI;
    let mu = mean - 0.5772 * beta;
    (mu, beta, mu - beta * (-(1.0 - ALPHA).ln()).ln())
}

/// Python's f"{x:.3e}" (two-digit signed exponent).
fn py_e3(x: f64) -> String {
    let s = format!("{:.3e}", x);
    let (m, e) = s.split_once('e').unwrap();
    let e: i32 = e.parse().unwrap();
    format!("{m}e{}{:02}", if e < 0 { '-' } else { '+' }, e.abs())
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.iter().any(|a| a == "--selftest") {
        std::process::exit(if nprng::selftest() { 0 } else { 1 });
    }
    let get = |k: &str, d: &str| -> String {
        args.iter().position(|a| a == k).map(|i| args[i + 1].clone()).unwrap_or_else(|| d.to_string())
    };
    let (fgp, bgp) = (&args[1], &args[2]);
    let perm: usize = get("--perm", "20").parse().unwrap();
    let seed: u64 = get("--seed", "271828").parse().unwrap();
    let threads: usize = get("--threads", "0").parse().unwrap();
    if threads > 0 {
        rayon::ThreadPoolBuilder::new().num_threads(threads).build_global().unwrap();
    }
    let t0 = Instant::now();
    let fg = read_fasta(fgp);
    let bg = read_fasta(bgp);
    let t_read = t0.elapsed().as_secs_f64();
    let shapes: Vec<Vec<Shape>> = (KMIN..=KMAX).map(all_shapes).collect();
    let lsh = shape(&(0..L).collect::<Vec<_>>());
    let dsh = shape(&[0, 1]);
    let t1 = Instant::now();
    let fr: Vec<&[u8]> = fg.iter().map(|v| v.as_slice()).collect();
    let br: Vec<&[u8]> = bg.iter().map(|v| v.as_slice()).collect();
    let (sf, sb) = (set_of(&fr), set_of(&br));
    let (cf, cb) = (count_all(&sf, &shapes, &lsh, &dsh), count_all(&sb, &shapes, &lsh, &dsh));
    let t_count = t1.elapsed().as_secs_f64();
    let t2 = Instant::now();
    let real = detect(&cf, &cb, &shapes, true);
    let t_stats = t2.elapsed().as_secs_f64();
    // permutation max-T: pool = fg + bg; only the permuted fg is counted
    let t3 = Instant::now();
    let pool: Vec<&[u8]> = fr.iter().chain(br.iter()).cloned().collect();
    let tot = Counts {
        g: cf.g.iter().zip(&cb.g).map(|(a, b)| a.iter().zip(b).map(|(u, v)| u.iter().zip(v).map(|(x, y)| x + y).collect()).collect()).collect(),
        l: cf.l.iter().zip(&cb.l).map(|(x, y)| x + y).collect(),
        di: cf.di.iter().zip(&cb.di).map(|(x, y)| x + y).collect(),
    };
    let mut rng = nprng::Pcg64::new(seed);      // numpy default_rng(271828), as proto/hlat.py
    let (mut mg, mut ml) = (Vec::new(), Vec::new());
    for _ in 0..perm {
        let idx = rng.permutation(pool.len());
        let p: Vec<&[u8]> = idx[..fg.len()].iter().map(|&i| pool[i]).collect();
        let cp = count_all(&set_of(&p), &shapes, &lsh, &dsh);
        let cq = sub(&tot, &cp);
        let r = detect(&cp, &cq, &shapes, false);
        if r.gapped_max_test.is_finite() {
            mg.push(r.gapped_max_test);
        }
        if r.lmer_max_test.is_finite() {
            ml.push(r.lmer_max_test);
        }
    }
    let t_perm = t3.elapsed().as_secs_f64();
    let fam_g: f64 = (KMIN..=KMAX).map(|k| all_shapes(k).len() as f64 * 4f64.powi(k as i32) / 2.0).sum();
    let fam_l = 4f64.powi(L as i32) / 2.0;
    let gg = if mg.len() >= 5 { Some(gumbel(&mg)) } else { None };
    let gl = if ml.len() >= 5 { Some(gumbel(&ml)) } else { None };
    let thr_g = gg.map_or(isf(0.025 / fam_g), |g| g.2);
    let thr_l = gl.map_or(isf(0.025 / fam_l), |g| g.2);
    let mut t_repr = 0.0;
    let out_prefix = get("-o", "");
    if !out_prefix.is_empty() {
        let t4 = Instant::now();
        let fgp_ = represent::padded(&fg);
        let bgp_ = represent::padded(&bg);
        let valid = |a: &represent::PArr| a.x.iter().filter(|&&c| c < 4).count() as f64;
        let r = valid(&fgp_) / valid(&bgp_);
        let (mf, mb) = (represent::markov1(&cf.di), represent::markov1(&cb.di));
        let sets = represent::make_sets(fgp_, bgp_, &mf, &mb);
        let seeds: Vec<represent::Seed> = real.seeds.iter().map(|(f, zr, zt, p)| represent::Seed {
            family: f.clone(), z_rank: *zr, z_test: *zt, pattern: p.clone() }).collect();
        let kept = represent::represent(&sets, &seeds, r, represent::GCRIT_001, 6.0);
        t_repr = t4.elapsed().as_secs_f64();
        if std::env::var("HASTUR_PROFILE").is_ok() {
            let tot: u64 = represent::PROF.iter().map(|a| a.load(std::sync::atomic::Ordering::Relaxed)).sum();
            eprintln!("representation wall {t_repr:.2}s; cumulative CPU by stage:");
            for (n, a) in represent::PROF_NAMES.iter().zip(represent::PROF.iter()) {
                let v = a.load(std::sync::atomic::Ordering::Relaxed) as f64 / 1e9;
                eprintln!("  {n:20} {v:8.2}s  {:5.1}%", 100.0 * v * 1e9 / tot.max(1) as f64);
            }
        }
        let mut tsv = String::from("rank\tconsensus\tz\tsig\tz_heldout\tbits\twidth\trepresentation\tsource\tseed\tp_adj\n");
        let mut meme = String::from("MEME version 4\n\nALPHABET= ACGT\n\nstrands: + -\n\nBackground letter frequencies\nA 0.25 C 0.25 G 0.25 T 0.25\n\n");
        for (i, d) in kept.iter().enumerate() {
            let (g, thr, fam) = if d.source == "gapped" { (gg, thr_g, fam_g) } else { (gl, thr_l, fam_l) };
            let sig = (d.z >= thr) as u8;
            let p_adj = match g {
                Some((mu, beta, _)) => (2.0 * -(-(-(d.z - mu) / beta.max(1e-9)).exp()).exp_m1()).min(1.0),
                None => (2.0 * fam * 0.5 * libm_erfc(d.z / std::f64::consts::SQRT_2)).min(1.0),
            };
            let iu = represent::iupac(&d.pwm);
            tsv += &format!("{}\t{}\t{:.2}\t{}\t{:.2}\t{:.1}\t{}\t{}\t{}\t{}\t{}\n", i + 1, iu, d.z, sig, d.z_heldout,
                            represent::info_bits(&d.pwm), d.pwm.len(), d.kind, d.source, d.seed, py_e3(p_adj));
            meme += &format!("MOTIF {}-{} hlat\nletter-probability matrix: alength= 4 w= {} nsites= {}\n", i + 1, iu, d.pwm.len(), d.nsites);
            for row in &d.pwm {
                meme += &format!("{:.4} {:.4} {:.4} {:.4}\n", row[0], row[1], row[2], row[3]);
            }
            meme += "\n";
        }
        std::fs::write(format!("{out_prefix}.tsv"), tsv).unwrap();
        std::fs::write(format!("{out_prefix}.meme"), meme).unwrap();
        std::fs::write(format!("{out_prefix}.json"), format!(
            "{{\"version\": \"v18-rust\", \"perm\": {perm}, \"n_fg\": {}, \"n_bg\": {}, \"zcrit\": {{\"gapped\": {thr_g}, \"lmer\": {thr_l}}}, \"perm_max\": {{\"gapped\": [{}], \"lmer\": [{}]}}, \"runtime_s\": {}}}",
            fg.len(), bg.len(), mg.iter().map(|x| x.to_string()).collect::<Vec<_>>().join(", "),
            ml.iter().map(|x| x.to_string()).collect::<Vec<_>>().join(", "), t0.elapsed().as_secs_f64())).unwrap();
    }
    let total = t0.elapsed().as_secs_f64();
    // JSON
    let seeds: Vec<String> = real.seeds.iter().map(|(f, zr, zt, s)| {
        let thr = if f == "gapped" { thr_g } else { thr_l };
        format!("{{\"family\":\"{f}\",\"pattern\":\"{s}\",\"z_rank\":{zr:.6},\"z_test\":{zt:.6},\"sig\":{}}}", (*zt >= thr) as u8)
    }).collect();
    let fmt = |v: &[f64]| v.iter().map(|x| format!("{x:.4}")).collect::<Vec<_>>().join(",");
    println!(
        "{{\"n_fg\":{},\"n_bg\":{},\"lambda_rank\":{{\"gapped\":{:.6},\"lmer\":{:.6}}},\"zcrit\":{{\"gapped\":{:.6},\"lmer\":{:.6}}},\
\"perm_max\":{{\"gapped\":[{}],\"lmer\":[{}]}},\"seeds\":[{}],\"time_s\":{{\"read\":{:.4},\"count\":{:.4},\"stats\":{:.4},\"perm\":{:.4},\"repr\":{:.4},\"total\":{:.4}}},\"threads\":{}}}",
        fg.len(), bg.len(), real.lam_rank_g, real.lam_rank_l, thr_g, thr_l, fmt(&mg), fmt(&ml), seeds.join(","),
        t_read, t_count, t_stats, t_perm, t_repr, total, rayon::current_num_threads()
    );
}
