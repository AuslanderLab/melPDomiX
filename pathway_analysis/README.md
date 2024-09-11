## Overview
This README describes the code used in Figure 3d-f and supplementary odds ratio plots.

## Calculation of pathway scores
**01_score_pathways.py**
Using a database of 31 pre-selected pathways, a pathway score for each patient is calculated for the omics types as follows:
- Whole-Exome Sequencing (WES): The number of genes in a pathway a patient has a mutation in are counted, and counts are zscored by pathway.
- RNA Sequencing (RNA-seq): For each patient and pathway, the zscored values for gene expression in a pathway are averaged.
- Reverse-Phase Protein Assay (RPPA): For each patient and pathway, the zscored values for protein abundance in a pathway are averaged.

**02_correlate_omics.py**
Pathway score dataframes for WES, RNA-seq, and RPPA are imported and labeled with their respective omics type. Clinical data is used to extract patient IDs for each phenotype (responder, nonresponder, mixed response/stable disease) for the respective treatment type (targeted therapy or immunotherapy). Individual omics dataframes are subsetted for patient IDs of each phenotype, and subsetted dataframes for each phenotype are concatenated across omics type. Correlation matrices are created for each phenotype using a Spearman correlation and finally saved in both 2D and 1D formats.

**03_plot_figures.R**
Pathways of interest are selected and reordered for ease of visualization. After slight dataframe reformatting, correlation matrices for each phenotype and treatment type are plotted using the corrplot package. Additional code is included for a strategy tested where RPPA data with isoforms of the same gene are first averaged before calculating overall RPPA score.

**Packages Used**
Pathway scores were calculated in Python 3.9.17 with the following packages:
- pandas
- scipy
Plotting was performed in R 4.2.1 using the following packages:
- dplyr
- ggplot2
- corrplot
