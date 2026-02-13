#PBS -q short
#PBS -l nodes=1:ppn=1
#PBS -l mem=8gb
#PBS -S /bin/bash
#PBS -N prinseq 
#PBS -j oe
#PBS -o log_prinseq

source /home/kn/kn_kn/kn_pop524719/.bashrc

mkdir -p trim2ends
cd trim2ends


prinseq++ -fastq ${libid}.sga4.filter.pass.fa.gz -out_name ${libid}.sga4.filter.pass.trim2ends.fq.gz \
 -min_len 30 -derep -trim_tail_left 6 -trim_tail_right 6
