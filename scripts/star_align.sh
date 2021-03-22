#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --time=72:00:00
#SBATCH --mem=120GB
#SBATCH --output=STAR_align.log
#SBATCH --job-name=STAR_align

module load intel-python3
source activate RNAseq_cp

job_start=`date +%s`

# Align RNA reads with STAR to the genome
# Memory and ulimit (# open files) errors occur if number of threads is too high
# Script can be improved to run all samples in parallel

for sample in $(cd ~/pw_oilfield/trimmed && ls *.fq | sed s/_[12]_trimmed.fq// | sort -u)
do

echo "Mapping ${sample}"

    STAR \
    --runThreadN 6 \
    --runMode alignReads \
    --genomeDir ~/pw_oilfield/assembly \
    --quantMode GeneCounts \
    --outSAMtype BAM SortedByCoordinate \
    --limitBAMsortRAM 32000000000 `#32GB - can probably be increased`\
    --readFilesIn ~/pw_oilfield/trimmed/${sample}_1_trimmed.fq ~/pw_oilfield/trimmed/${sample}_2_trimmed.fq \
    --outFileNamePrefix ~/pw_oilfield/alignment_sorted/${sample}

done

echo "Finished STAR align"
job_end=`date +%s`
runtime=$((job_end-job_start))
hours=$((runtime / 3600))
minutes=$(( (runtime % 3600) / 60 ))
seconds=$(( (runtime % 3600) % 60 ))
echo "seqkit stats runtime: $hours:$minutes:$seconds (hh:mm:ss)"
