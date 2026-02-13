# shotgun data processing
This repository contains scripts for intial data processing of shotgun sequenced DNA libraries.  
Step 1: Use [fastp](https://github.com/OpenGene/fastp) to remove adapters, read filtering by length and get statistics of sequence libraries.  
Step 2: use [SGA](https://github.com/jts/sga) to remove exact duplicates and low complexity reads.  
Step 3: Use [Prinseq++](https://github.com/Adrian-Cantu/PRINSEQ-plus-plus) to trim Poly A/T head and tails.  
