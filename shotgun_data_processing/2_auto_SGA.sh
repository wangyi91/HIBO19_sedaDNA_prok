#PBS -q tiny
#PBS -l nodes=1:ppn=1
#PBS -l mem=1gb
#PBS -S /bin/bash
#PBS -N auto
#PBS -j oe
#PBS -o LOG_auto_SGA

source ~/.bashrc
cd $hb/fastp

for file in $hb/fastp/*.trimmed.merged.fastq.gz
do 
loc=${file%.trimmed.merged.fastq.gz}
libid=${loc#*fastp/}
libnum=${libid#AVXF-1-}
qsub -v loc=$loc,libid=$libid -N sga_$libnum -o log_sga4_$libnum $hb/sga.sh
done
