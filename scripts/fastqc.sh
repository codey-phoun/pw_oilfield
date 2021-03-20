#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=28
#SBATCH --time=72:00:00
#SBATCH --mem-per-cpu=4500
#SBATCH --output=fastqc.log
#SBATCH --job-name=fastqc

module load intel-python3
source activate RNAseq_v2.0

job_start=`date +%s`

fastqc ~/pw_oilfield/trimmed/*.fq.gz \
--threads 56 \
--outdir ~/pw_oilfield/fastqc

multiqc ~/pw_oilfield/fastqc ~/pw_oilfield/trimmed \
--filename fastqc_multiqc_trimmed.html \
--outdir ~/pw_oilfield/fastqc/multiqc_trimmed

echo "Finished FastQC"
job_end=`date +%s`
runtime=$((job_end-job_start))
hours=$((runtime / 3600))
minutes=$(( (runtime % 3600) / 60 ))
seconds=$(( (runtime % 3600) % 60 ))
echo "FastQC runtime: $hours:$minutes:$seconds (hh:mm:ss)"