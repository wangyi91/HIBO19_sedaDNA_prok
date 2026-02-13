#PBS -q short
#PBS -l nodes=1:ppn=14
#PBS -l mem=32gb
#PBS -S /bin/bash
#PBS -N fastp
#PBS -j oe


source /home/kn/kn_kn/kn_pop524719/.bashrc
mkdir -p $hb/fastp
cd $hb/fastp

conda activate base

:<<'END'

# option 1: output merged fastq
fastp -w 14 --in1 ${loc}_R1.fastq.gz --in2 ${loc}_R2.fastq.gz --merge --merged_out $hb/fastp/${libid}.trimmed.merged.fastq.gz \
 --out1 $hb/fastp/${libid}_R1.trimmed.unmerged.fastq.gz --out2 $hb/fastp/${libid}_R2.trimmed.unmerged.fastq.gz \
 -V --detect_adapter_for_pe --dedup --dup_calc_accuracy 6 \
 -e 25 -y -p -x --poly_x_min_len 7 \
 -q 30 -l 30 -h $hb/fastp/${libid}.html

END

# option 2: output filtered read1 and read2 separately
fastp -w 14 --in1 ${loc}_R1.fastq.gz --in2 ${loc}_R2.fastq.gz \
--out1 $hb/fastp/${libid}_R1.trimmed.paired.fastq.gz --out2 $hb/fastp/${libid}_R2.trimmed.paired.fastq.gz \
 --unpaired1 $hb/fastp/${libid}_R1.trimmed.unpaired.fastq.gz --unpaired2 $hb/fastp/${libid}_R2.trimmed.unpaired.fastq.gz \
 -V --detect_adapter_for_pe --dedup --dup_calc_accuracy 6 \
 -e 25 -y -p -x --poly_x_min_len 7 \
 -q 30 -l 30 -h $hb/fastp/${libid}.html

# option 3 following 2: trim 3bp at each end, and output filtered read1 and read2 separately





:<<'END'

# bulk processing all libraries
for infile in $data/*${N}_R1.fastq.gz
do
loc=${infile%_R1*}
libid=${loc#*L002_}
AdapterRemoval --file1 "${loc}_R1.fastq.gz"  --file2 "${loc}_R2.fastq.gz" --mm 3 --minlength 24 \
 --basename $libid --trimns --trimqualities --minquality 30 \
 --adapter1 AATGATACGGCGACCACCGAGATCTACACNNNNNNNACACTCTTTCCCTACACGACGCTCTTCCGATCT \
 --adapter2 CAAGCAGAAGACGGCATACGAGATNNNNNNNGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT \
 --collapse
done


# merge collapsed files
for file in *collapsed
do cat $file ${file}.truncated > ${file}.merged.fq
done

END
