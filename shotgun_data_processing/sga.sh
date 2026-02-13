#PBS -q short
#PBS -l nodes=1:ppn=14
#PBS -l mem=8gb
#PBS -S /bin/bash
#PBS -N sga
#PBS -j oe

source ~/.bashrc
cd $hb/fastp

sga preprocess -m 30 --dust-threshold=4 ${loc}.trimmed.merged.fastq.gz -o ${loc}.sga4.fastq
sga index --algorithm=ropebwt --threads=14 ${loc}.sga4.fastq
sga filter --verbose --threads=14 --low-complexity-check --no-kmer-check ${loc}.sga4.fastq

wc -l ${libid}.sga4.discard.fa > lines.${libid}.sga4.discard
wc -l ${libid}.sga4.filter.pass.fa > lines.${libid}.sga4.filter.pass


