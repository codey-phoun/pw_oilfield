#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --time=72:00:00
#SBATCH --mem=8GB
#SBATCH --output=STAR_align_parallel_submit.log
#SBATCH --job-name=STAR_align_submit

module load intel-python3
source activate RNAseq_cp

# Submits multiple STAR jobs with SLURM to align each sample to the genome

for sample in $(cd ~/pw_oilfield/trimmed && ls *.fq | sed s/_[12]_trimmed.fq// | sort -u)
do

sbatch ~/pw_oilfield/scripts/star_align_parallel.sh ${sample}
echo "Submitted ${sample}"

done

echo "Finished STAR align parallel submission"
