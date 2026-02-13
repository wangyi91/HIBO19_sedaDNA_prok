#PBS -q tiny
#PBS -l nodes=1:ppn=1
#PBS -l mem=8gb
#PBS -S /bin/bash
#PBS -N auto
#PBS -j oe
#PBS -o LOG_auto_fastp

source /home/kn/kn_kn/kn_pop524719/.bashrc
mkdir -p $hb/fastp
cd $hb/fastp

for file in $data/*_R1.fastq.gz
do
loc=${file%_R1*}
libid=${loc#*L002_}
libnum=${libid#AVXF-1-}

qsub -v loc=$loc,libid=$libid -N fastp_$libnum -o log_fastp_$libnum $hb/scripts/fastp.sh
done
