# Produced Water - Oilfield : RNA Seq Analysis of *Phaeodactylum tricornutum*

By:

Codey Phoun

Claudia Vo

Sudhin Domala

&nbsp;

For:

SJSU CS286 Project

Dr. Andreopoulos 

Spring 2021

&nbsp;

## Project Introduction

Produced water (PW) is water produced as a byproduct created during the extraction of oil and natural gases. The water quality of PW ranges from well to well, but most PW contains oils, heavy metals, and traces of naturally occurring radioactive material. PW can serve as an alternative water source if the pollutants found in the water are removed. Bioremediation by algae is one potential method for the treatment of  PW. Alternatively, the PW may potentially be used as a growth medium for algae cultivation in biofuel production.

*Phaeodactylum tricornutum* is a model marine diatom with a fully sequenced genome. This diatom was previously observed by Dr. Jeroen Gillard to have increased photosynthetic capacity when grown in pure PW compared to a normal medium. An RNA sequencing project was conducted by his lab to investigate the differential gene expression of *Phaeodactylum tricornutum* grown under three different growth medium conditions and different growth phases. In this study, we processed the RNA sequencing data and explored which Phaeodactylum tricornutum genes were differentially expressed under the different growth mediums and growth phases.

&nbsp;

## Project Files

RNA_seq_project file contains the data processing steps of the RNA seq data. All shell scripts used for processing on the SJSU CoS HPC are located in the scripts/ directory. The STAR_results/ directory contains the gene counts matrix created by merging the STAR results. The R/ contains the R markdown file, html and pdf render of the markdown file, and the output .xlsx files of the DEGs and gene ontology results. The html render of the R markdown file can also be viewed at: <https://codey-phoun.github.io/pw_oilfield/>
