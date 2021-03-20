#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=28
#SBATCH --time=72:00:00
#SBATCH --mem-per-cpu=4500
#SBATCH --output=cutadapt.log
#SBATCH --job-name=cutadapt

module load intel-python3
source activate RNAseq_cp

start=`date +%s`

for sample in HQ_ST1 HQ_ST2 PW_ST1 PW_ST2 HQE1 HQE2 HQ10E1 HQ10E2 PWE1 PWE2

do
    echo "Processing $sample"

    cutadapt \
    -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA `# TruSeq single index Read 1` \
    -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT `# # TruSeq single index Read 2` \
    --cores 0 `# auto-detect number of CPU cores to use` \
    --overlap 10 `# minimum length overlap between read and adapter` \
    --quality-cutoff 20,20 `# trim low-quality bases from 5' and 3' ends before adapter removal` \
    --minimum-length 50 `# minimum read length after trimming` \
    -o ~/pw_oilfield/trimmed/${sample}_1_trimmed.fq.gz `# Read 1 output` \
    -p ~/pw_oilfield/trimmed/${sample}_2_trimmed.fq.gz `# Read 2 output` \
    ~/pw_oilfield/data/${sample}_1.fq.gz `# Read 1 input` \
    ~/pw_oilfield/data/${sample}_2.fq.gz `# Read 2 input` \
    > ~/pw_oilfield/trimmed/${sample}_trim_log.txt

done
echo "Finished adapter trimming"

# Run FastQC and MultiQC on trimmed reads
fastqc ~/pw_oilfield/trimmed/*.fq.gz \
--threads 56 \
--outdir ~/pw_oilfield/fastqc
echo "Finished FastQC"

multiqc ~/pw_oilfield/fastqc ~/pw_oilfield/trimmed \
--filename fastqc_multiqc_trimmed.html \
--outdir ~/pw_oilfield/fastqc/multiqc_trimmed
echo "Finished MultiQC"

end=`date +%s`
runtime=$((end-start))
hours=$((runtime / 3600))
minutes=$(( (runtime % 3600) / 60 ))
seconds=$(( (runtime % 3600) % 60 ))
echo "Runtime: $hours:$minutes:$seconds (hh:mm:ss)"