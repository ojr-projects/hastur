# %% run_tools.py: run every tool on every dataset directory (fg.fa + bg.fa), cached, with timeouts
"""
Tools (same fg and bg for all):
  hlat        this project's prototype (proto/hlat.py); hlat_f2r = frozen F2 from the Rust binary (HASTUR_BIN)
  meme        MEME ZOOPS, order-1 background from bg (E1); xstreme = XSTREME, no motif database (E1)
  streme      STREME 5.5.9, --nmotifs 10
  homer       HOMER findMotifs.pl fasta mode, -len 8,10,12
  comse       Comsequent v1.0 all-in-one (discover k5/w13 + extend -join 30 with theta tables), as in its P5
Outputs go to <dataset>/<tool>/; a tool is skipped when its output exists. Timeouts kill the whole process group.
Usage: python bench/run_tools.py <root> [tools=hlat,streme,homer,comse] [parallel=48]
"""
import os
import signal
import subprocess
import sys
import time
from multiprocessing import Pool

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
MM = os.path.expanduser('~/micromamba/envs')
PY = os.path.join(MM, 'bench', 'bin', 'python')
COMSE = os.path.expanduser(os.environ.get('COMSE_BIN', '~/comsequent/comse'))
TIMEOUT = int(os.environ.get('TOOL_TIMEOUT', '7200'))
HASTUR_BIN = os.path.expanduser(os.environ.get('HASTUR_BIN', '~/.cache/hastur_target/release/hastur'))


MIN_FREE_GB = float(os.environ.get('MIN_FREE_GB', '4'))
LAST_RC = {}


def mem_available_gb():
    for line in open('/proc/meminfo'):
        if line.startswith('MemAvailable'):
            return int(line.split()[1]) / 1048576
    return 0.0


def timed(cmd, log, env=None):
    while mem_available_gb() < MIN_FREE_GB:      # memory gate: a 15 GB workstation OOM-killed STREME (2026-10-10)
        time.sleep(15)
    t0 = time.time()
    with open(log, 'w') as fh:
        p = subprocess.Popen(cmd, stdout=fh, stderr=subprocess.STDOUT, start_new_session=True, env=env)
        try:
            p.wait(timeout=TIMEOUT)
        except subprocess.TimeoutExpired:
            os.killpg(p.pid, signal.SIGKILL)
            p.wait()
            fh.write(f'\nTIMEOUT after {TIMEOUT} s\n')
    LAST_RC[log] = p.returncode
    return time.time() - t0


def run_one(args):
    d, tool = args
    fg, bg = os.path.join(d, 'fg.fa'), os.path.join(d, 'bg.fa')
    o = os.path.join(d, tool)
    os.makedirs(o, exist_ok=True)
    done = os.path.join(o, 'runtime_s')
    if os.path.exists(done):
        return f'{os.path.basename(d)} {tool} cached'
    env = dict(os.environ, NUMBA_NUM_THREADS='2', OMP_NUM_THREADS='1', PATH=os.path.join(MM, tool.split('_')[0], 'bin') + ':' + os.environ['PATH'])
    if tool.startswith('hlat') and tool != 'hlat_f2r':        # HASTUR arms: hlatc = --center 40 (LN1); hlat_a1/a2/a12 = LN2 variants
        extra = {'hlat': [], 'hlatc': ['--center', '40'], 'hlat_a1': ['--sites'], 'hlat_a2': ['--gtrim', '0.01'],
                 'hlat_a12': ['--sites', '--gtrim', '0.01'], 'hlat_a3': ['--sites', '--gtrim', '0.01'],
                 'hlat_c1': ['--repr3', '--gtrim', '0.01'], 'hlat_c2': ['--sites', '--gtrim', '0.01', '--seeds', '20'],
                 'hlat_c12': ['--repr3', '--gtrim', '0.01', '--seeds', '20'],
                 'hlat_e1': ['--sites', '--gtrim', '0.01', '--exact'],
                 'hlat_l1': ['--sites', '--gtrim', '0.01', '--refit', '2'],                   # LN5
                 'hlat_l2': ['--sites', '--gtrim', '0.01', '--plainz'],
                 'hlat_l12': ['--sites', '--gtrim', '0.01', '--refit', '2', '--plainz'],
                 'hlat_x1': ['--sites', '--gtrim', '0.01', '--exact-all'],                     # LN6
                 'hlat_x2': ['--sites', '--gtrim', '0.01', '--lmer-valid'],
                 'hlat_x12': ['--sites', '--gtrim', '0.01', '--exact-all', '--lmer-valid'],
                 'hlat_y1': ['--sites', '--gtrim', '0.01', '--exact-all', '--seqcount'],          # LN6 round 2
                 'hlat_y2': ['--sites', '--gtrim', '0.01', '--exact-all', '--seqcount', '--lmer-ball'],
                 'hlat_y3': ['--sites', '--gtrim', '0.01', '--exact-all', '--lmer-ball'],
                 'hlat_z1': ['--sites', '--gtrim', '0.01', '--exact-all', '--seqcount', '--lmer-ball', '--presence-emp'],   # LN6 r3
                 'hlat_p2': ['--sites', '--gtrim', '0.01', '--exact-all', '--lmer-ball', '--perm', '20'],          # LN7 dev
                 'hlat_p3': ['--sites', '--gtrim', '0.01', '--exact-all', '--seqcount', '--lmer-ball', '--perm', '20'],
                 'hlat_p5': ['--sites', '--gtrim', '0.01', '--exact-all', '--lmer-ball', '--perm', '20', '--rank-v8'],
                 'hlat_p6': ['--sites', '--gtrim', '0.01', '--exact-all', '--lmer-ball', '--perm', '20', '--rank-v8', '--perm-bg'],
                 'hlat_f1': ['--sites', '--gtrim', '0.01', '--plainz', '--exact-all', '--lmer-ball', '--perm', '20', '--rank-v8'],
                 'hlat_p7': ['--sites', '--gtrim', '0.01', '--exact-all', '--lmer-ball', '--perm', '20', '--rank-v8', '--no-gc-test'],
                 'hlat_f2': ['--sites', '--gtrim', '0.01', '--plainz', '--exact-all', '--lmer-ball', '--perm', '20', '--rank-v8', '--no-gc-test'],
                 'hlat_u1': ['--sites', '--gtrim', '0.01', '--plainz', '--exact-all', '--lmer-ball', '--perm', '20', '--rank-v8', '--no-gc-test', '--union']}[tool]
        dt = timed([PY, os.path.join(REPO, 'proto', 'hlat.py'), '-i', fg, '-b', bg, '-o', os.path.join(o, 'hlat')] + extra,
                   os.path.join(o, 'log'), env)
    elif tool == 'hlat_f2r':        # frozen F2 from the Rust binary (byte-identical to hlat_f2; ROADMAP_v2 rule 6)
        dt = timed([HASTUR_BIN, fg, bg, '-o', os.path.join(o, 'hlat'), '--threads', os.environ.get('HASTUR_THREADS', '2')],
                   os.path.join(o, 'log'), env)
    elif tool == 'meme':            # MEME ZOOPS (E1): order-1 background from bg, both strands, 5 motifs reported
        mk = os.path.join(o, 'bg.markov')
        dt = timed([os.path.join(MM, 'meme', 'bin', 'fasta-get-markov'), '-m', '1', '-dna', bg, mk], os.path.join(o, 'mklog'), env)
        dt += timed([os.path.join(MM, 'meme', 'bin', 'meme'), fg, '-dna', '-oc', o, '-mod', 'zoops', '-nmotifs', '5',
                     '-minw', '6', '-maxw', '25', '-revcomp', '-bfile', mk, '-nostatus'], os.path.join(o, 'log'), env)
    elif tool == 'xstreme':         # XSTREME (E1): STREME + MEME with the same background, no motif database
        dt = timed([os.path.join(MM, 'meme', 'bin', 'xstreme'), '--p', fg, '--n', bg, '--dna', '--oc', o],
                   os.path.join(o, 'log'), env)
    elif tool == 'sea':             # SEA known-motif enrichment with the bank in SEA_BANK (E1 known-motif arm)
        dt = timed([os.path.join(MM, 'meme', 'bin', 'sea'), '--p', fg, '--n', bg, '--m', os.environ['SEA_BANK'], '--text'],
                   os.path.join(o, 'sea.tsv'), env)    # --text: the table on stdout (this install lacks the HTML template)
    elif tool == 'streme':
        dt = timed([os.path.join(MM, 'meme', 'bin', 'streme'), '--p', fg, '--n', bg, '--dna', '--nmotifs', '10', '--oc', o],
                   os.path.join(o, 'log'), env)
    elif tool == 'streme_all':      # FR1: report everything, thresholds swept afterwards
        dt = timed([os.path.join(MM, 'meme', 'bin', 'streme'), '--p', fg, '--n', bg, '--dna', '--thresh', '1',
                    '--nmotifs', '25', '--oc', o], os.path.join(o, 'log'), env)
    elif tool == 'homer':
        dt = timed([os.path.join(MM, 'homer', 'bin', 'findMotifs.pl'), fg, 'fasta', o, '-fasta', bg, '-noknown', '-p', '2',
                    '-len', '8,10,12'], os.path.join(o, 'log'), env)
    elif tool == 'comse':
        pre = os.path.join(o, 'cs')
        dt = timed([COMSE, 'discover', '-i', fg, '-bg', bg, '-o', pre, '-k', '5', '-w', '13', '-strand', 'both', '-t', '2',
                    '-arrays=false'], pre + '.log')
        dt += timed([COMSE, 'extend', '-i', fg, '-bg', bg, '-seeds', pre + '.calls_enr_tmax.tsv', '-top', '10', '-kmax', '14',
                     '-wmax', '20', '-iupac', '-join', '30', '-t', '2', '-theta', pre + '.calls_theta.tsv', '-thetamax',
                     pre + '.calls_theta_tmax.tsv', '-o', pre + '.aio.tsv'], pre + '.aio.log')
    else:
        raise ValueError(tool)
    rcs = [rc for lg, rc in LAST_RC.items() if lg.startswith(o)]
    if tool == 'hlat_f2r' and not os.path.exists(os.path.join(o, 'hlat.tsv')):    # a stale binary ran detection only
        rcs.append('no hlat.tsv')
    if any(rc != 0 for rc in rcs):      # killed (e.g. OOM) or failed: never leave a completion marker
        with open(os.path.join(o, 'failed'), 'w') as fh:
            fh.write(f'return codes {rcs}\n')
        return f'{os.path.basename(d)} {tool} FAILED rc={rcs}'
    with open(done, 'w') as fh:
        fh.write(f'{dt:.1f}\n')
    return f'{os.path.basename(d)} {tool} {dt:.0f}s'


if __name__ == '__main__':
    root = sys.argv[1]
    tools = (sys.argv[2] if len(sys.argv) > 2 else 'hlat,streme,homer,comse').split(',')
    par = int(sys.argv[3]) if len(sys.argv) > 3 else 48
    ds = sorted(os.path.join(root, x) for x in os.listdir(root) if os.path.exists(os.path.join(root, x, 'fg.fa')))
    jobs = [(d, t) for t in tools for d in ds]
    with Pool(par) as pool:
        for msg in pool.imap_unordered(run_one, jobs):
            print(msg, flush=True)
