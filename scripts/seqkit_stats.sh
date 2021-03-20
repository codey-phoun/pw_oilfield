#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=28
#SBATCH --time=72:00:00
#SBATCH --mem-per-cpu=4500
#SBATCH --output=seqkit.log
#SBATCH --job-name=seqkit

module load intel-python3
source activate RNAseq_cp

job_start=`date +%s`
seqkit stats -Ta -j 28 ~/pw_oilfield/data/*.fq.gz | csvtk csv2md -t > ~/pw_oilfield/data/data_file_stats.txt

echo "Finished seqkit stats"
job_end=`date +%s`
runtime=$((job_end-job_start))
hours=$((runtime / 3600))
minutes=$(( (runtime % 3600) / 60 ))
seconds=$(( (runtime % 3600) % 60 ))
echo "seqkit stats runtime: $hours:$minutes:$seconds (hh:mm:ss)"