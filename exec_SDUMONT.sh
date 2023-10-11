#!/bin/bash

module load python/3.8.2
echo 'Python LOADED'

declare -A datasets=(
    ["./B.1.427_429/B.1.427_B.1.429-USA-California_MSA.afa"]="B.1.427_429"
    ["./B.1.525/B.1.525-UK-Nigeria_MSA.afa"]="B.1.525"
    ["./B.1.526/BB.1.526-NewYork_MSA.afa"]="B.1.526"
    ["./B.1.617.1/B.1.617.1-India_MSA.afa"]="B.1.617.1"
    ["./C.37/C.37-Peru_MSA.afa"]="C.37"
    ["./P.1/P.1-Brazil-Japan_MSA.afa"]="P.1"
    ["./P.2/P.2-Brazil_MSA.afa"]="P.2"
    ["./P.3/P.3-Philippines_MSA.afa"]="P.3"
)

# Loop through the datasets
for path in "${!datasets[@]}"; do
    base="${datasets[$path]}"
    srun --nodes=1 --ntasks=1 --ntasks-per-node=1 --cpus-per-task=24 -p dockthor python3 consenso.py "$path" "./$base/${base}_cosensus.csv" "./$base/${base}_canonical_cosensus.fasta"
    echo "${base} FINISH"
done

# Special case for B.1.1.7
for i in {1..6}; do
    srun --nodes=1 --ntasks=1 --ntasks-per-node=1 --cpus-per-task=24 -p dockthor python3 consenso.py "./B.1.1.7/$i/B.1.1.7-UK-pt${i}_MSA.afa" "./B.1.1.7/$i/B.1.1.7-UK-pt${i}_cosensus.csv" "./B.1.1.7/$i/B.1.1.7-UK-pt${i}_canonical_cosensus.fasta"
    echo "B.1.1.7-UK-pt${i} FINISH"
done
