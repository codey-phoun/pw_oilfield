#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=28
#SBATCH --time=72:00:00
#SBATCH --mem-per-cpu=4500
#SBATCH --output=cutadapt_test.log
#SBATCH --job-name=cutadapt

module load intel-python3
source activate RNAseq_cp

job_start=`date +%s`

cutadapt \
-b AGATCGGAAGAGCACACGTCTGAACTCCAGTCA `# LTHT_1` \
-b TGACTGGAGTTCAGACGTGTGCTCTTCCGATCT `# LTHT_1_RC` \
-b AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT `# LTHT_2` \
-b ACACTCTTTCCCTACACGACGCTCTTCCGATCT `# LTHT_2_RC` \
-b TGGAATTCTCGGGTGCCAAGG `# SmallRNA` \
-b CCTTGGCACCCGAGAATTCCA `# SmallRNA_RC` \
-b CTGTCTCTTATACACATCT `# Nextera_1` \
-b AGATGTGTATAAGAGACAG `# Nextera_2/Nextera_1_RC`\
-b AATGATACGGCGACCACCGAGATCT `# TruSeq Universal Adapter` \
-b ATCTCGTATGCCGTCTTCTGCTTG `# TruSeq Index Adapter` \
-B AGATCGGAAGAGCACACGTCTGAACTCCAGTCA `# LTHT_1` \
-B TGACTGGAGTTCAGACGTGTGCTCTTCCGATCT `# LTHT_1_RC` \
-B AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT `# LTHT_2` \
-B ACACTCTTTCCCTACACGACGCTCTTCCGATCT `# LTHT_2_RC` \
-B TGGAATTCTCGGGTGCCAAGG `# SmallRNA` \
-B CCTTGGCACCCGAGAATTCCA `# SmallRNA_RC` \
-B CTGTCTCTTATACACATCT `# Nextera_1` \
-B AGATGTGTATAAGAGACAG `# Nextera_2/Nextera_1_RC`\
-B AATGATACGGCGACCACCGAGATCT `# TruSeq Universal Adapter` \
-B ATCTCGTATGCCGTCTTCTGCTTG `# TruSeq Index Adapter` \
--cores 0 `# auto-detect number of CPU cores to use` \
--overlap 15 `# minimum length overlap between read and adapter` \
--times 10 `# max number of adapters removed from each read` \
--quality-cutoff 20,20 `# trim low-quality bases from 5' and 3' ends before adapter removal` \
-o ~/pw_oilfield/trimmed/HQ_ST1_1_trimmed.fq `# Read 1 output` \
-p ~/pw_oilfield/trimmed/HQ_ST1_2_trimmed.fq `# Read 2 output` \
~/pw_oilfield/data/HQ_ST1_1.fq `# Read 1 input` \
~/pw_oilfield/data/HQ_ST1_2.fq `# Read 2 input` \
> ~/pw_oilfield/trimmed/HQ_ST1_trim_log.txt

echo "Finished adapter trimming"
end=`date +%s`
runtime=$((end-start))
hours=$((runtime / 3600))
minutes=$(( (runtime % 3600) / 60 ))
seconds=$(( (runtime % 3600) % 60 ))
echo "Runtime: $hours:$minutes:$seconds (hh:mm:ss)"