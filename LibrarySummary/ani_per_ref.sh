#!/usr/bin/env bash

BAM_DIR="../../EnvironmenTracker/results_ani70/taxonomic-profiling"
OUT_DIR="./input"
mkdir -p "$OUT_DIR"

for bam in "${BAM_DIR}"/HB_*.dedup.bam; do
    lib=$(basename "$bam" .dedup.bam)

    echo "Processing $lib ..."

    samtools view -F 0x904 "$bam" | \
    awk '{
        seq_len = length($10)
        nm = 0
        for(i=12;i<=NF;i++) {
            if($i ~ /^NM:i:/) {
                split($i,a,":"); nm=a[3]
            }
        }
        pid = (seq_len - nm) / seq_len * 100
        ref = $3
        sum[ref] += pid
        count[ref]++
    }
    END {
        for(ref in sum) {
            ani = sum[ref] / count[ref]
            print ref, ani, count[ref]
        }
    }' > "${OUT_DIR}/${lib}.ANI_per_ref.txt"

    echo "${lib}.ANI_per_ref.txt written to ${OUT_DIR}/"
done

