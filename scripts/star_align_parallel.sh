#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --time=72:00:00
#SBATCH --mem=120GB
#SBATCH --job-name=STAR_align

module load intel-python3
source activate RNAseq_cp

job_start=`date +%s`

# Align RNA reads with STAR to the genome
# Memory and ulimit (# open files) errors occur if number of threads is too high

sample=$1

echo "Mapping ${sample}"

STAR \
--runThreadN 12 \
--runMode alignReads \
--genomeDir ~/pw_oilfield/assembly \
--quantMode GeneCounts \
--outSAMtype BAM SortedByCoordinate \
--limitBAMsortRAM 32000000000 `#32GB - can be increased if needed`\
--readFilesIn ~/pw_oilfield/trimmed/${sample}_1_trimmed.fq ~/pw_oilfield/trimmed/${sample}_2_trimmed.fq \
--outFileNamePrefix ~/pw_oilfield/alignment_sorted/${sample}_

echo "Finished STAR align"
job_end=`date +%s`
runtime=$((job_end-job_start))
hours=$((runtime / 3600))
minutes=$(( (runtime % 3600) / 60 ))
seconds=$(( (runtime % 3600) % 60 ))
echo "STAR align runtime: $hours:$minutes:$seconds (hh:mm:ss)"
