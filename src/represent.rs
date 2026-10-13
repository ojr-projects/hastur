//! Representation stage of HASTUR F2 (proto/hlat.py run(): widened frame, held-out z of frame vs EM, EM ZOOPS,
//! trimming, keep / de-duplicate, site refit with the plain threshold), ported for identical output.

use crate::nprng::Pcg64;
use rayon::prelude::*;
use std::sync::atomic::{AtomicU64, Ordering};
use std::time::Instant;

/// Cumulative CPU time per stage (ns), for profiling: 0 widened_frame, 1 frame_init (held-out fits), 2 EM window
/// scoring, 3 EM posterior/exp/accumulate, 4 scores() (held-out + site), 5 site best windows, 6 other.
pub static PROF: [AtomicU64; 7] = [AtomicU64::new(0), AtomicU64::new(0), AtomicU64::new(0), AtomicU64::new(0),
                                     AtomicU64::new(0), AtomicU64::new(0), AtomicU64::new(0)];
#[inline]
fn prof(i: usize, t: Instant) {
    PROF[i].fetch_add(t.elapsed().as_nanos() as u64, Ordering::Relaxed);
}
pub const PROF_NAMES: [&str; 7] = ["widened_frame", "frame_init_heldout", "em_scoring", "em_posterior", "scores", "site_best_windows", "other"];

const LN_1E6_SAFE: f64 = -13.815510557964274 - 1e-9; // ln(1e-6) minus a margin for exp's rounding
pub const GCRIT_001: f64 = 11.344866730144368; // scipy.stats.chi2.isf(0.01, 3)
const G_CRIT: f64 = 16.27;

/// numpy's pairwise summation (float add.reduce over a contiguous 1-D array).
pub fn np_sum(a: &[f64]) -> f64 {
    let n = a.len();
    if n < 8 {
        let mut r = 0.0;
        for x in a {
            r += x;
        }
        r
    } else if n <= 128 {
        let mut r = [a[0], a[1], a[2], a[3], a[4], a[5], a[6], a[7]];
        let mut i = 8;
        while i < n - n % 8 {
            for j in 0..8 {
                r[j] += a[i + j];
            }
            i += 8;
        }
        let mut res = ((r[0] + r[1]) + (r[2] + r[3])) + ((r[4] + r[5]) + (r[6] + r[7]));
        while i < n {
            res += a[i];
            i += 1;
        }
        res
    } else {
        let mut n2 = n / 2;
        n2 -= n2 % 8;
        np_sum(&a[..n2]) + np_sum(&a[n2..])
    }
}

// ---------------------------------------------------------------- padded arrays (Python's 2-D uint8 matrices)
#[derive(Clone)]
pub struct PArr {
    pub n: usize,
    pub l: usize,
    pub x: Vec<u8>,
}

impl PArr {
    #[inline]
    pub fn row(&self, i: usize) -> &[u8] {
        &self.x[i * self.l..(i + 1) * self.l]
    }
}

pub fn padded(seqs: &[Vec<u8>]) -> PArr {
    let l = seqs.iter().map(|s| s.len()).max().unwrap_or(0);
    let mut x = vec![4u8; seqs.len() * l];
    for (i, s) in seqs.iter().enumerate() {
        x[i * l..i * l + s.len()].copy_from_slice(s);
    }
    PArr { n: seqs.len(), l, x }
}

pub fn revcomp(a: &PArr) -> PArr {
    let mut x = vec![4u8; a.x.len()];
    for i in 0..a.n {
        let r = a.row(i);
        let mut last = a.l;
        while last > 0 && r[last - 1] == 4 {
            last -= 1;
        }
        for j in 0..last {
            let c = r[last - 1 - j];
            x[i * a.l + j] = if c < 4 { 3 - c } else { 4 };
        }
    }
    PArr { n: a.n, l: a.l, x }
}

#[derive(Clone)]
pub struct Model1 {
    pub pi: [f64; 4],
    pub t: [[f64; 4]; 4],
}

/// Python markov(X, R, 1): (2-mer counts both strands) + 0.5; pi = row sums / total (numpy pairwise total).
pub fn markov1(di: &[f64]) -> Model1 {
    let c: Vec<f64> = di.iter().map(|v| v + 0.5).collect();
    let tot = np_sum(&c);
    let mut pi = [0.0; 4];
    let mut t = [[0.0; 4]; 4];
    for a in 0..4 {
        let rs = ((c[a * 4] + c[a * 4 + 1]) + c[a * 4 + 2]) + c[a * 4 + 3];
        pi[a] = rs / tot;
        for b in 0..4 {
            t[a][b] = c[a * 4 + b] / rs;
        }
    }
    Model1 { pi, t }
}

/// Python bg_logprob for order 1: log pi (renormalised) at t = 0, log T[x_{t-1}, x_t] after; N/padding count as T.
pub fn bg_logprob(a: &PArr, m: &Model1) -> Vec<f64> {
    let s = ((m.pi[0] + m.pi[1]) + m.pi[2]) + m.pi[3];
    let lp0: Vec<f64> = m.pi.iter().map(|p| (p / s).ln()).collect();
    let mut lt = [[0.0; 4]; 4];
    for i in 0..4 {
        for j in 0..4 {
            lt[i][j] = m.t[i][j].ln();
        }
    }
    let mut out = vec![0.0; a.x.len()];
    for i in 0..a.n {
        let r = a.row(i);
        if a.l == 0 {
            continue;
        }
        out[i * a.l] = lp0[r[0].min(3) as usize];
        for t in 1..a.l {
            out[i * a.l + t] = lt[r[t - 1].min(3) as usize][r[t].min(3) as usize];
        }
    }
    out
}

/// Python sample_chain (order 1) with numpy's Generator draws.
pub fn sample_chain(a: &PArr, m: &Model1, rng: &mut Pcg64) -> PArr {
    let (n, l) = (a.n, a.l);
    let mut y = vec![0u8; n * l];
    let d = rng.choice_p(&m.pi, n);
    let mut st: Vec<usize> = d.clone();
    for i in 0..n {
        if l > 0 {
            y[i * l] = d[i] as u8;
        }
    }
    let mut cdf = [[0.0; 4]; 4];
    for s in 0..4 {
        let mut acc = 0.0;
        for b in 0..4 {
            acc += m.t[s][b];
            cdf[s][b] = acc;
        }
    }
    for t in 1..l {
        let u: Vec<f64> = (0..n).map(|_| rng.random()).collect();
        for i in 0..n {
            let v = cdf[st[i]].iter().filter(|&&c| u[i] > c).count();
            y[i * l + t] = v as u8;
            st[i] = (st[i] * 4 + v) % 4;
        }
    }
    for k in 0..y.len() {
        if a.x[k] > 3 {
            y[k] = 4;
        }
    }
    PArr { n, l, x: y }
}

/// One set: sequences, reverse complements, per-position background log-probabilities of both under its own chain.
pub struct RSet {
    pub x: PArr,
    pub r: PArr,
    pub g0: Vec<f64>,
    pub g1: Vec<f64>,
    pub len: Vec<usize>, // row length without trailing padding (same on both strands)
}

pub fn rset(x: PArr, m: &Model1) -> RSet {
    let r = revcomp(&x);
    let g0 = bg_logprob(&x, m);
    let g1 = bg_logprob(&r, m);
    let len = (0..x.n).map(|i| {
        let row = x.row(i);
        let mut last = x.l;
        while last > 0 && row[last - 1] == 4 {
            last -= 1;
        }
        last
    }).collect();
    RSet { x, r, g0, g1, len }
}

/// Python _llr_windows: per-window log-likelihood ratios of sequence i, both strands, -inf where a window holds N.
/// Column-outer / window-inner, so the inner loop vectorises; each window still accumulates its columns in Python's
/// order (v = 0; v += lp[j][c] - G[t+j] for j = 0..W), so every value is bit-identical. An N maps to a -inf table
/// entry, which leaves the window at exactly -inf, as Python's early exit does.
fn llr_windows(s: &RSet, i: usize, lp: &[[f64; 4]], out: &mut Vec<f64>) {
    let w = lp.len();
    let l = s.x.l;
    out.clear();
    if l < w {
        return;
    }
    let lw = l - w + 1;
    // windows starting at t > len - w reach the trailing padding: exactly -inf in Python (early exit on code 4)
    let nv = if s.len[i] >= w { s.len[i] - w + 1 } else { 0 };
    out.resize(2 * lw, 0.0);
    for strand in 0..2 {
        let (y, g) = if strand == 0 { (&s.x, &s.g0) } else { (&s.r, &s.g1) };
        let row = y.row(i);
        let gr = &g[i * l..(i + 1) * l];
        let o = &mut out[strand * lw..(strand + 1) * lw];
        for j in 0..w {
            let (t0, t1, t2, t3) = (lp[j][0], lp[j][1], lp[j][2], lp[j][3]);
            let rj = &row[j..j + nv];
            let gj = &gr[j..j + nv];
            let tab = [t0, t1, t2, t3, f64::NEG_INFINITY];
            for t in 0..nv {
                o[t] += tab[rj[t].min(4) as usize] - gj[t];
            }
        }
        for t in nv..lw {
            o[t] = f64::NEG_INFINITY;
        }
    }
}

/// Best window of sequence i: (score, strand, position); numpy argmax (first maximum; -inf if all -inf).
fn best_window(s: &RSet, i: usize, lp: &[[f64; 4]], buf: &mut Vec<f64>) -> (f64, usize, usize) {
    llr_windows(s, i, lp, buf);
    if buf.is_empty() {
        return (f64::NEG_INFINITY, 0, 0);
    }
    let mut bi = 0;
    for k in 1..buf.len() {
        if buf[k] > buf[bi] {
            bi = k;
        }
    }
    let lw = buf.len() / 2;
    if bi < lw { (buf[bi], 0, bi) } else { (buf[bi], 1, bi - lw) }
}

pub fn scores(s: &RSet, idx: &[usize], lp: &[[f64; 4]]) -> Vec<f64> {
    let t = Instant::now();
    let r = idx.par_iter().map_init(Vec::new, |buf, &i| best_window(s, i, lp, buf).0).collect();
    prof(4, t);
    r
}

// ---------------------------------------------------------------- matrix helpers
pub type Pwm = Vec<[f64; 4]>;

fn rowsum(r: &[f64; 4]) -> f64 {
    ((r[0] + r[1]) + r[2]) + r[3]
}

pub fn column_g(p: &Pwm, n: f64) -> Vec<f64> {
    p.iter().map(|r| 2.0 * n * (((r[0] * (r[0] / 0.25).ln() + r[1] * (r[1] / 0.25).ln()) + r[2] * (r[2] / 0.25).ln()) + r[3] * (r[3] / 0.25).ln())).collect()
}

pub fn trim_sig(p: &Pwm, n: f64, gcrit: f64) -> Pwm {
    let g = column_g(p, n);
    let keep: Vec<usize> = (0..g.len()).filter(|&i| g[i] >= gcrit).collect();
    if keep.is_empty() { p[..1].to_vec() } else { p[keep[0]..=*keep.last().unwrap()].to_vec() }
}

fn ic_rows(p: &Pwm) -> Vec<f64> {
    p.iter().map(|r| 2.0 + (((r[0] * r[0].log2() + r[1] * r[1].log2()) + r[2] * r[2].log2()) + r[3] * r[3].log2())).collect()
}

pub fn trim(p: &Pwm) -> Pwm {
    let ic = ic_rows(p);
    let keep: Vec<usize> = (0..ic.len()).filter(|&i| ic[i] >= 0.3).collect();
    if keep.is_empty() { p[..1].to_vec() } else { p[keep[0]..=*keep.last().unwrap()].to_vec() }
}

pub fn info_bits(p: &Pwm) -> f64 {
    np_sum(&ic_rows(p))
}

pub fn pwm_pwm_sim(a: &Pwm, b: &Pwm) -> f64 {
    let (a, b) = if a.len() > b.len() { (b, a) } else { (a, b) };
    let u = [0.25; 4];
    let brev: Pwm = b.iter().rev().map(|r| [r[3], r[2], r[1], r[0]]).collect();
    let mut best = 0.0f64;
    for m in [b, &brev] {
        let (la, lm) = (a.len() as isize, m.len() as isize);
        for off in -la + 1..lm {
            let mut s = 0.0;
            for j in 0..a.len() {
                let i = j as isize + off;
                let col = if i >= 0 && i < lm { &m[i as usize] } else { &u };
                let tv = (((col[0] - a[j][0]).abs() + (col[1] - a[j][1]).abs()) + (col[2] - a[j][2]).abs()) + (col[3] - a[j][3]).abs();
                s += 1.0 - 0.5 * tv;
            }
            best = best.max(s / a.len() as f64);
        }
    }
    best
}

pub fn iupac(p: &Pwm) -> String {
    let mut out = String::new();
    for col in p {
        let mut o = [0usize, 1, 2, 3];
        o.sort_by(|&a, &b| col[b].partial_cmp(&col[a]).unwrap()); // stable: ties keep index order (numpy n=4)
        let ch = if col[o[0]] >= 0.6 {
            b"ACGT"[o[0]] as char
        } else if col[o[0]] + col[o[1]] >= 0.8 {
            let (x, y) = (o[0].min(o[1]), o[0].max(o[1]));
            match (x, y) { (0, 1) => 'M', (0, 2) => 'R', (0, 3) => 'W', (1, 2) => 'S', (1, 3) => 'Y', _ => 'K' }
        } else {
            'N'
        };
        out.push(ch);
    }
    let t = out.trim_matches('N');
    if t.is_empty() { "N".to_string() } else { t.to_string() }
}

// ---------------------------------------------------------------- numpy quantile grid and contrasts
pub fn quantile_grid(vals: &[f64]) -> Vec<f64> {
    let mut v: Vec<f64> = vals.iter().cloned().filter(|x| x.is_finite()).collect();
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let n = v.len();
    let step = (0.98 - 0.3) / 24.0;
    (0..25).map(|i| {
        let q = if i == 24 { 0.98 } else { i as f64 * step + 0.3 };
        if n == 0 {
            return f64::NAN;
        }
        let nf = n as f64;
        let vi = nf * q + (1.0 + q * (1.0 - 1.0 - 1.0)) - 1.0; // numpy _compute_virtual_index, alpha = beta = 1
        let prev = vi.floor();
        let gamma = vi - prev;
        let pi = (prev as isize).clamp(0, n as isize - 1) as usize;
        let ni = (prev as isize + 1).clamp(0, n as isize - 1) as usize;
        let (a, b) = (v[pi], v[ni]);
        let d = b - a;
        if gamma >= 0.5 { b - d * (1.0 - gamma) } else { a + d * gamma }
    }).collect()
}

/// dc_z over the four sets in Python's order (f, b, sf, sb).
pub fn dc_z(sc: &[Vec<f64>; 4], t: f64) -> f64 {
    let mut a = [0.0; 4];
    let mut rate = [0.0; 4];
    for k in 0..4 {
        let c = sc[k].iter().filter(|&&v| v >= t).count() as f64 + 0.5;
        let n = sc[k].len() as f64 + 1.0;
        a[k] = c;
        rate[k] = c / n;
    }
    let d = (rate[0] / rate[2]).ln() - (rate[1] / rate[3]).ln();
    let mut s = 0.0;
    for k in 0..4 {
        s += (1.0 - rate[k]) / a[k];
    }
    d / s.sqrt()
}

pub fn plain_z(sc: &[Vec<f64>; 4], t: f64) -> f64 {
    let a = sc[0].iter().filter(|&&v| v >= t).count() as f64;
    let na = sc[0].len() as f64;
    let b = sc[1].iter().filter(|&&v| v >= t).count() as f64;
    let nb = sc[1].len() as f64;
    let p = (a + b) / (na + nb);
    (a / na - b / nb) / ((p * (1.0 - p)).max(1e-12) * (1.0 / na + 1.0 / nb)).sqrt()
}

fn argmax_first<F: Fn(f64) -> f64>(grid: &[f64], f: F) -> f64 {
    let mut best = grid[0];
    let mut bv = f(grid[0]);
    for &g in &grid[1..] {
        let v = f(g);
        if v > bv {
            bv = v;
            best = g;
        }
    }
    best
}

// ---------------------------------------------------------------- frame and EM
pub fn frame_init(f: &RSet, fi: &[usize], b: &RSet, bi: &[usize], offs: &[usize], letters: &[u8], flank: isize, r: f64) -> (Pwm, f64) {
    let span = offs[offs.len() - 1] + 1;
    let lf = f.x.l as isize;
    let fl = 0.max(flank.min((lf - span as isize).div_euclid(2))) as usize;
    let w = span + 2 * fl;
    let mut cols = Vec::new();
    for (s, idx) in [(f, fi), (b, bi)] {
        let mut cnt = vec![[0.0f64; 4]; w];
        for y in [&s.x, &s.r] {
            if y.l < w {
                continue;
            }
            let per: Vec<Vec<[f64; 4]>> = idx.par_iter().map(|&i| {
                let row = y.row(i);
                let mut c = vec![[0.0f64; 4]; w];
                for p in 0..=(y.l - w) {
                    if offs.iter().zip(letters).all(|(&o, &lt)| row[p + fl + o] == lt) {
                        for j in 0..w {
                            let v = row[p + j];
                            if v < 4 {
                                c[j][v as usize] += 1.0;
                            }
                        }
                    }
                }
                c
            }).collect();
            for c in per {
                for j in 0..w {
                    for k in 0..4 {
                        cnt[j][k] += c[j][k];
                    }
                }
            }
        }
        cols.push(cnt);
    }
    let raw: Vec<[f64; 4]> = (0..w).map(|j| {
        let mut r4 = [0.0; 4];
        for k in 0..4 {
            r4[k] = cols[0][j][k] - r * cols[1][j][k];
        }
        r4
    }).collect();
    let m: Pwm = raw.iter().map(|r4| {
        let ex = [r4[0].max(0.0) + 0.5, r4[1].max(0.0) + 0.5, r4[2].max(0.0) + 0.5, r4[3].max(0.0) + 0.5];
        let s = rowsum(&ex);
        [ex[0] / s, ex[1] / s, ex[2] / s, ex[3] / s]
    }).collect();
    let rs: Vec<f64> = raw.iter().map(rowsum).collect();
    let n = (np_sum(&rs) / w as f64).max(1.0);
    (m, n)
}

pub fn widened_frame(f: &RSet, fi: &[usize], b: &RSet, bi: &[usize], offs: &[usize], letters: &[u8], r: f64) -> (Pwm, isize, f64) {
    let span = offs[offs.len() - 1] + 1;
    let mut flank: isize = 4;
    loop {
        let (m, n) = frame_init(f, fi, b, bi, offs, letters, flank, r);
        let g = column_g(&m, n);
        let at_edge = g[0] >= G_CRIT || g[g.len() - 1] >= G_CRIT;
        if !at_edge || span as isize + 2 * (flank + 4) > 30.min(f.x.l as isize) {
            return (m, flank, n);
        }
        flank += 4;
    }
}

pub fn em_zoops(s: &RSet, idx: &[usize], pwm0: &Pwm, gamma0: f64, iters: usize, pseudo: f64) -> (Pwm, f64) {
    let mut pwm = pwm0.clone();
    let mut gamma = gamma0;
    let w = pwm.len();
    let l = s.x.l;
    let n = idx.len() as f64;
    let mut sc: Vec<f64> = Vec::new();
    for _ in 0..iters {
        let lp: Vec<[f64; 4]> = pwm.iter().map(|r| [r[0].ln(), r[1].ln(), r[2].ln(), r[3].ln()]).collect();
        let mut cnt = vec![[0.0f64; 4]; w];
        let mut tot = 0.0;
        for &i in idx {
            let t1 = Instant::now();
            llr_windows(s, i, &lp, &mut sc);
            prof(2, t1);
            let t2 = Instant::now();
            let m = sc.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
            if !m.is_finite() {
                continue;
            }
            let mut nw = 0usize;
            let mut acc = 0.0;
            for &v in &sc {
                if v.is_finite() {
                    nw += 1;
                    acc += (v - m).exp();
                }
            }
            let lse = m + acc.ln();
            let a1 = (1.0 - gamma).ln();
            let a2 = (gamma / nw as f64).ln() + lse;
            let logden = a1.max(a2) + (-(a1 - a2).abs()).exp().ln_1p();
            let lw = l - w + 1;
            let lgn = (gamma / nw as f64).ln();
            for (k, &v) in sc.iter().enumerate() {
                if !v.is_finite() {
                    continue;
                }
                let e = lgn + v - logden;
                if e < LN_1E6_SAFE {
                    continue; // exp(e) < 1e-6 for certain: Python computes it and discards the window
                }
                let pst = e.exp();
                if pst < 1e-6 {
                    continue;
                }
                tot += pst;
                let (y, t) = if k < lw { (&s.x, k) } else { (&s.r, k - lw) };
                let row = y.row(i);
                for j in 0..w {
                    cnt[j][row[t + j] as usize] += pst;
                }
            }
            prof(3, t2);
        }
        let new: Pwm = cnt.iter().map(|c| {
            let rs = ((c[0] + c[1]) + c[2]) + c[3];
            [(c[0] + pseudo * 0.25) / (rs + pseudo), (c[1] + pseudo * 0.25) / (rs + pseudo),
             (c[2] + pseudo * 0.25) / (rs + pseudo), (c[3] + pseudo * 0.25) / (rs + pseudo)]
        }).collect();
        gamma = (tot / n).max(0.02).min(0.98);
        let mut diff = 0.0f64;
        for j in 0..w {
            for k in 0..4 {
                diff = diff.max((new[j][k] - pwm[j][k]).abs());
            }
        }
        pwm = new;
        if diff < 1e-4 {
            break;
        }
    }
    (pwm, gamma)
}

// ---------------------------------------------------------------- the four sets and held-out evidence
pub struct Sets {
    pub s: [RSet; 4], // f, b, sf, sb
    pub all: [Vec<usize>; 4],
    pub half: [[Vec<usize>; 2]; 4],  // rows with i % 2 == h (test)
    pub train: [[Vec<usize>; 2]; 4], // rows with i % 2 != h
}

pub fn make_sets(fg: PArr, bg: PArr, mf: &Model1, mb: &Model1) -> Sets {
    let mut rng = Pcg64::new(0);
    let sf = sample_chain(&fg, mf, &mut rng);
    let sb = sample_chain(&bg, mb, &mut rng);
    let s = [rset(fg, mf), rset(bg, mb), rset(sf, mf), rset(sb, mb)];
    let all = [0, 1, 2, 3].map(|k| (0..s[k].x.n).collect::<Vec<usize>>());
    let half = [0, 1, 2, 3].map(|k| [0, 1].map(|h| (0..s[k].x.n).filter(|i| i % 2 == h).collect::<Vec<usize>>()));
    let train = [0, 1, 2, 3].map(|k| [0, 1].map(|h| (0..s[k].x.n).filter(|i| i % 2 != h).collect::<Vec<usize>>()));
    Sets { s, all, half, train }
}

fn sc4(sets: &Sets, idx: [&Vec<usize>; 4], lp: &[[f64; 4]]) -> [Vec<f64>; 4] {
    [0, 1, 2, 3].map(|k| scores(&sets.s[k], idx[k], lp))
}

pub fn heldout_z<F: Fn(usize) -> Pwm>(sets: &Sets, fit: F) -> f64 {
    let mut z = Vec::new();
    for h in 0..2 {
        let p = trim(&fit(h));
        let lp: Vec<[f64; 4]> = p.iter().map(|r| [r[0].ln(), r[1].ln(), r[2].ln(), r[3].ln()]).collect();
        let tr = [&sets.train[0][h], &sets.train[1][h], &sets.train[2][h], &sets.train[3][h]];
        let te = [&sets.half[0][h], &sets.half[1][h], &sets.half[2][h], &sets.half[3][h]];
        let sc_tr = sc4(sets, tr, &lp);
        let grid = quantile_grid(&sc_tr[0]);
        let t = argmax_first(&grid, |t| dc_z(&sc_tr, t));
        z.push(dc_z(&sc4(sets, te, &lp), t));
    }
    (0.0 + z[0] + z[1]) / std::f64::consts::SQRT_2
}

/// Site refit (F2: plain threshold, all data, one round). Returns (matrix, nsites).
pub fn site_matrix(sets: &Sets, m: &Pwm, gcrit: f64) -> (Pwm, usize) {
    let lp: Vec<[f64; 4]> = m.iter().map(|r| [r[0].ln(), r[1].ln(), r[2].ln(), r[3].ln()]).collect();
    let all = [&sets.all[0], &sets.all[1], &sets.all[2], &sets.all[3]];
    let sc = sc4(sets, all, &lp);
    let grid = quantile_grid(&sc[0]);
    let t = argmax_first(&grid, |t| plain_z(&sc, t));
    let f = &sets.s[0];
    let t5 = Instant::now();
    let bw: Vec<(f64, usize, usize)> = (0..f.x.n).into_par_iter().map_init(Vec::new, |buf, i| best_window(f, i, &lp, buf)).collect();
    prof(5, t5);
    let sel: Vec<usize> = (0..bw.len()).filter(|&i| bw[i].0 >= t).collect();
    if sel.len() < 10 {
        return (m.clone(), sel.len());
    }
    let w = m.len();
    let mut cnt = vec![[0.0f64; 4]; w];
    for &i in &sel {
        let (_, st, ps) = bw[i];
        let row = if st == 0 { f.x.row(i) } else { f.r.row(i) };
        for j in 0..w {
            let c = row[ps + j];
            if c < 4 {
                cnt[j][c as usize] += 1.0;
            }
        }
    }
    let out: Pwm = cnt.iter().map(|c| {
        let rs = rowsum(c) + 1.0;
        [(c[0] + 0.25) / rs, (c[1] + 0.25) / rs, (c[2] + 0.25) / rs, (c[3] + 0.25) / rs]
    }).collect();
    let out = trim_sig(&out, sel.len() as f64, gcrit);
    (if pwm_pwm_sim(&out, m) >= 0.6 { out } else { m.clone() }, sel.len())
}

// ---------------------------------------------------------------- the per-seed pipeline and the keep loop
pub struct Seed {
    pub family: String,
    pub z_rank: f64,
    pub z_test: f64,
    pub pattern: String,
}

pub struct Kept {
    pub z: f64,
    pub z_heldout: f64,
    pub kind: String,
    pub source: String,
    pub seed: String,
    pub pwm: Pwm,
    pub nsites: usize,
}

struct SeedFit {
    pwm: Pwm,
    kind: &'static str,
    z_heldout: f64,
    n_sites: f64,
}

fn fit_seed(sets: &Sets, sd: &Seed, r: f64, gcrit: f64) -> SeedFit {
    let offs: Vec<usize> = sd.pattern.bytes().enumerate().filter(|(_, c)| *c != b'N').map(|(i, _)| i).collect();
    let letters: Vec<u8> = sd.pattern.bytes().filter(|c| *c != b'N').map(|c| match c { b'A' => 0, b'C' => 1, b'G' => 2, _ => 3 }).collect();
    let (f, b) = (&sets.s[0], &sets.s[1]);
    let t0 = Instant::now();
    let (init, fl, n_sites) = widened_frame(f, &sets.all[0], b, &sets.all[1], &offs, &letters, r);
    prof(0, t0);
    let fit_frame = |h: usize| { let t = Instant::now(); let m = frame_init(f, &sets.train[0][h], b, &sets.train[1][h], &offs, &letters, fl, r).0; prof(1, t); m };
    let fit_em = |h: usize| em_zoops(f, &sets.train[0][h], &fit_frame(h), 0.3, 60, 0.25).0;
    let z_frame = heldout_z(sets, fit_frame);
    let z_em = heldout_z(sets, fit_em);
    let mut kind = if z_em > z_frame { "em" } else { "frame" };
    let mut pwm = trim_sig(&init, n_sites, gcrit);
    if kind == "em" {
        let (em_m, em_g) = em_zoops(f, &sets.all[0], &init, 0.3, 60, 0.25);
        let em_all = trim_sig(&em_m, em_g * f.x.n as f64, gcrit);
        if pwm_pwm_sim(&em_all, &pwm) >= 0.6 {
            pwm = em_all;
        } else {
            kind = "frame";
        }
    }
    SeedFit { pwm, kind, z_heldout: z_em.max(z_frame), n_sites }
}

pub fn represent(sets: &Sets, seeds: &[Seed], r: f64, gcrit: f64, min_bits: f64) -> Vec<Kept> {
    // Python processes seeds in order of -z_rank (stable); each seed's fit is independent of the others, so all
    // fits run in parallel and the order-dependent keep / de-duplicate pass runs afterwards, as in Python.
    let mut order: Vec<usize> = (0..seeds.len()).collect();
    order.sort_by(|&a, &b| seeds[b].z_rank.partial_cmp(&seeds[a].z_rank).unwrap());
    let fits: Vec<SeedFit> = order.par_iter().map(|&i| fit_seed(sets, &seeds[i], r, gcrit)).collect();
    let mut kept: Vec<(Kept, Pwm)> = Vec::new();
    for (k, &i) in order.iter().enumerate() {
        let sd = &seeds[i];
        // explaining away: consistent with an already-kept motif's selection-stage consensus
        if kept.iter().any(|(_, sel)| {
            let iu: String = iupac(sel).chars().map(|c| if "ACGT".contains(c) { c } else { 'N' }).collect();
            crate::consistent_overlap(&sd.pattern, &iu, 4)
        }) {
            continue;
        }
        let ft = &fits[k];
        if info_bits(&ft.pwm) < min_bits || kept.iter().any(|(_, sel)| pwm_pwm_sim(&ft.pwm, sel) >= 0.75) {
            continue;
        }
        let pwm_sel = ft.pwm.clone();
        let (pwm, nsites) = site_matrix(sets, &ft.pwm, gcrit);
        let _ = ft.n_sites;
        kept.push((Kept { z: sd.z_test, z_heldout: ft.z_heldout, kind: format!("{}+sites", ft.kind), source: sd.family.clone(),
                          seed: sd.pattern.clone(), pwm, nsites }, pwm_sel));
    }
    kept.into_iter().map(|(k, _)| k).collect()
}
