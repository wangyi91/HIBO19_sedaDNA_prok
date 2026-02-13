#PBS -q tiny
#PBS -l nodes=1:ppn=1
#PBS -l mem=8gb
#PBS -S /bin/bash
#PBS -N auto
#PBS -j oe
#PBS -o LOG_auto_prinseq

source /home/kn/kn_kn/kn_pop524719/.bashrc
mkdir -p $hb/trim2ends
cd $hb/trim2ends

for file in $hb/sga/*.trimmed.merged.fastq.gz
do
loc=${file%_R1*}
libid=${loc#*L002_}
libnum=${libid#AVXF-1-}

qsub -v loc=$loc,libid=$libid -N prinseq_$libnum -o log_prinseq_$libnum $hb/scripts/prinseq.sh
done

