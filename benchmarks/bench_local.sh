#!/bin/bash
# Speed floor: Rust detection core vs Python F2 on the same machine. Python timings: total run and detection share
# (cumulative find_seeds over the real pass + 20 permutations, from cProfile).
P=$HOME/micromamba/envs/bench/bin/python
RB=$HOME/.cache/hastur_target/release/hastur
T=/tmp/claude-1000/-home-oliver-code-idea-coins/e99b72c4-978c-4d0d-a83a-9b1a8ca38f17/scratchpad/bench
mkdir -p $T
F2="--sites --gtrim 0.01 --plainz --exact-all --lmer-ball --perm 20 --rank-v8 --no-gc-test"
printf "%-24s %6s %6s | %9s %9s %9s | %8s %8s\n" unit n_fg n_bg py_total py_detect py_det_24t rust_1t rust_24t
for u in $HOME/runs/upnull_dev3/saureus_n20__r1 $HOME/runs/regulon/CRP $HOME/runs/lowN/CTCF__n250__r0 $HOME/runs/lowN/CTCF__n1000__r0 $HOME/data/chip/CTCF; do
    b=$(basename $u); nf=$(grep -c '>' $u/fg.fa); nb=$(grep -c '>' $u/bg.fa)
    # warm numba cache once
    NUMBA_NUM_THREADS=2 $P proto/hlat.py -i $u/fg.fa -b $u/bg.fa -o $T/warm $F2 > /dev/null 2>&1
    pyt=() 
    for th in 2 24; do
        s=$(date +%s.%N)
        NUMBA_NUM_THREADS=$th $P -c "
import cProfile, pstats, sys; sys.path.insert(0, 'proto'); import hlat
cProfile.run(\"hlat.run('$u/fg.fa', '$u/bg.fa', '$T/o', sites=True, gtrim_p=0.01, plainz=True, exact_all=True, lmer_ball=True, perm=20, rank_v8=True, no_gc_test=True, verbose=False)\", '$T/p')
st = pstats.Stats('$T/p'); print([v[3] for k, v in st.stats.items() if k[2] == 'find_seeds'][0])" > $T/det_$th.txt 2>/dev/null
        e=$(echo "$(date +%s.%N) - $s" | bc); pyt+=("$e" "$(cat $T/det_$th.txt)")
    done
    r1=$($RB $u/fg.fa $u/bg.fa --perm 20 --threads 1 | $P -c "import json,sys; print(json.load(sys.stdin)['time_s']['total'])")
    r24=$($RB $u/fg.fa $u/bg.fa --perm 20 --threads 24 | $P -c "import json,sys; print(json.load(sys.stdin)['time_s']['total'])")
    printf "%-24s %6s %6s | %8.1fs %8.1fs %8.1fs | %7.2fs %7.2fs\n" $b $nf $nb ${pyt[0]} ${pyt[1]} ${pyt[3]} $r1 $r24
done
