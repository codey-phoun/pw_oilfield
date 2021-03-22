# modified from https://ucdavis-bioinformatics-training.github.io/2020-mRNA_Seq_Workshop/data_reduction/03-counts_mm

# create header file
echo gene_name $(cd ~/pw_oilfield/alignment_sorted && ls *_ReadsPerGene.out.tab | sed s/_ReadsPerGene.out.tab// | sort -u) > ~/pw_oilfield/alignment_sorted/tmp/header.txt

# Place each sample's STAR gene count file - ReadsPerGene.out.tab in the tmp/ directory 
# The 2nd column (-f2) of ReadsPerGene.out.tab contains the non-stranded counts
for sample in $(cd ~/pw_oilfield/alignment_sorted && ls *_ReadsPerGene.out.tab | sed s/_ReadsPerGene.out.tab// | sort -u)
do 
    echo ${sample}
    cat ~/pw_oilfield/alignment_sorted/${sample}_ReadsPerGene.out.tab | tail -n +5 | cut -f2 > ~/pw_oilfield/alignment_sorted/tmp/${sample}.count
done

# get a list of gene ids (-f1)
tail -n +5 ~/pw_oilfield/alignment_sorted/HQ10E1_ReadsPerGene.out.tab | cut -f1 > ~/pw_oilfield/alignment_sorted/tmp/geneids.txt

# combine all the columns of the count files
paste ~/pw_oilfield/alignment_sorted/tmp/geneids.txt ~/pw_oilfield/alignment_sorted/tmp/*.count > ~/pw_oilfield/alignment_sorted/tmp/tmp.out

# add the header
cat <(cat ~/pw_oilfield/alignment_sorted/tmp/header.txt | sed 's/ /\t/g') ~/pw_oilfield/alignment_sorted/tmp/tmp.out > ~/pw_oilfield/alignment_sorted/STAR_counts.txt

# remove the tmp folder
rm -rf ~/pw_oilfield/alignment_sorted/tmp
