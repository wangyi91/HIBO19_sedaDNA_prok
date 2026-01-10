#!/usr/bin/bash
for bam in ../../EnvironmenTracker/results/taxonomic-profiling/HB_*.filtered.bam
do
    lib=$(basename "$bam" .dedup.filtered.bam)
    samtools view -F 0x904 "$bam" | awk '!seen[$1]++ { print length($10) }' > "./input/${lib}.filtered_bam_read_len.txt"
    echo ${lib}.filtered_bam_read_len.txt written to ./input/
done
