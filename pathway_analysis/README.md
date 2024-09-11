## Overview
This README describes the code used in Figure 3d-f and supplementary odds ratio plots.

## Calculation of pathway scores
**01_score_pathways.py**
Using a database of 31 pre-selected pathways, a pathway score for each patient is calculated for the omics types as follows:
- Whole-Exome Sequencing (WES): The number of genes in a pathway a patient has a mutation in are counted, and counts are zscored by pathway.
- RNA Sequencing (RNA-seq): For each patient and pathway, the zscored values for gene expression in a pathway are averaged.
- Reverse-Phase Protein Assay (RPPA): For each patient and pathway, the zscored values for protein abundance in a pathway are averaged.

**02_correlate_omics.py**


Pathway scores were calculated in Python 3.9.17 with the following packages:
- pandas
- scipy
- 
