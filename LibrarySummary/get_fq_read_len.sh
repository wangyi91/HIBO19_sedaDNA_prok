#!/usr/bin/bash
for fq in ../../trim2ends/*fq.gz
do
    lib="HB_$(basename "$fq" .sga4.filter.pass.trim2ends.fq.gz | sed 's/.*-//')"
    # Extract every 4th line starting from the second, then print sequence lengths
    zcat "$fq" | awk 'NR % 4 == 2 { print length($0) }' > "./input/${lib}.all_read_len.txt"

    echo "${lib}.all_read_len.txt written to ./input/"
done

