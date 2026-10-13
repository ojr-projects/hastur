#!/bin/bash
# Full pipeline speed (same outputs, validated): Python F2 vs Rust hastur, idle workstation.
P=$HOME/micromamba/envs/bench/bin/python
RB=$HOME/.cache/hastur_target/release/hastur
T=/tmp/claude-1000/-home-oliver-code-idea-coins/e99b72c4-978c-4d0d-a83a-9b1a8ca38f17/scratchpad/bfull
mkdir -p $T
F2="--sites --gtrim 0.01 --plainz --exact-all --lmer-ball --perm 20 --rank-v8 --no-gc-test"
printf "%-22s %6s %6s | %10s %10s | %9s %9s | %s\n" unit n_fg n_bg py_2thr py_24thr rust_1t rust_24t "rust 24t breakdown (count/perm/repr)"
for u in $HOME/runs/upnull_dev3/saureus_n20__r1 $HOME/runs/regulon/CRP $HOME/runs/lowN/CTCF__n250__r0 $HOME/runs/lowN/CTCF__n1000__r0 $HOME/data/chip/CTCF; do
    b=$(basename $u); nf=$(grep -c '>' $u/fg.fa); nb=$(grep -c '>' $u/bg.fa)
    tp=()
    for th in 2 24; do s=$(date +%s.%N); NUMBA_NUM_THREADS=$th $P proto/hlat.py -i $u/fg.fa -b $u/bg.fa -o $T/py $F2 > /dev/null 2>&1; tp+=($(echo "$(date +%s.%N) - $s" | bc)); done
    tr=()
    for th in 1 24; do s=$(date +%s.%N); $RB $u/fg.fa $u/bg.fa -o $T/rs --perm 20 --threads $th > $T/rs_$th.json; tr+=($(echo "$(date +%s.%N) - $s" | bc)); done
    bd=$($P -c "import json; t=json.load(open('$T/rs_24.json'))['time_s']; print(f\"{t['count']:.2f}/{t['perm']:.2f}/{t['repr']:.2f}\")")
    printf "%-22s %6s %6s | %9.1fs %9.1fs | %8.2fs %8.2fs | %s\n" $b $nf $nb ${tp[0]} ${tp[1]} ${tr[0]} ${tr[1]} $bd
done
