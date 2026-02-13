#PBS -q short
#PBS -l nodes=1:ppn=14
#PBS -l mem=8gb
#PBS -S /bin/bash
#PBS -N sga
#PBS -j oe

source /home/kn/kn_kn/kn_pop524719/.bashrc
cd $hb/fastp
conda activate shotgun

#loc='/beegfs/work/kn_pop524719/HIBO_shotgun/fastp/AVXF-1-30'
#libid='AVXF-1-30'

sga preprocess -m 30 --dust-threshold=4 ${loc}.trimmed.merged.fastq.gz -o ${loc}.sga4.fastq
sga index --algorithm=ropebwt --threads=14 ${loc}.sga4.fastq
sga filter --verbose --threads=14 --low-complexity-check --no-kmer-check ${loc}.sga4.fastq

wc -l ${libid}.sga4.discard.fa > lines.${libid}.sga4.discard
wc -l ${libid}.sga4.filter.pass.fa > lines.${libid}.sga4.filter.pass


