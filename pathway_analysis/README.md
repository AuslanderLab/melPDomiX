## Overview
This README describes the code used in Figure 3d-f and supplementary odds ratio plots.

## Calculation of Pathway Scores
**01_score_pathways.py**
Using a database of 31 pre-selected pathways, a pathway score for each patient is calculated for the omics types as follows:
- Whole-Exome Sequencing (WES): The number of genes in a pathway a patient has a mutation in are counted, and counts are zscored by pathway.
- RNA Sequencing (RNA-seq): For each patient and pathway, the zscored values for gene expression in a pathway are averaged.
- Reverse-Phase Protein Assay (RPPA): For each patient and pathway, the zscored values for protein abundance in a pathway are averaged. Only the first listed isoform of a gene is included in the average protein abundance for a pathway.

**02_correlate_omics.py**
Pathway score dataframes for WES, RNA-seq, and RPPA are imported and labeled with their respective omics type. Clinical data is used to extract patient IDs for each phenotype (responder, nonresponder, or mixed response/stable disease) for the respective treatment type (targeted therapy or immunotherapy). Individual omics dataframes are subsetted for patient IDs of each phenotype, and subsetted dataframes for each phenotype are concatenated across omics type. Correlation matrices are created for each phenotype using a Spearman correlation and finally saved in both 2D and 1D formats.

**03_plot_figures.R**
Pathways of interest are selected and reordered for ease of visualization. After slight dataframe reformatting, correlation matrices for each phenotype and treatment type are plotted using the corrplot package. Additional code is included for a strategy tested where RPPA data with isoforms of the same gene are first averaged before calculating overall RPPA score. Code for plotting the following correlation plots is included:
Patients given immunotherapy:
**- those who did not respond to treatment (figure 3d)**
**- those with mixed response to treatment or who maintained stable disease (figure 3e)**
**- those responded partially or completely to treatment (figure 3f)**

Patients given targeted therapy:
- those who did not respond to treatment
- those with mixed response to treatment or who maintained stable disease
- those responded partially or completely to treatment

**Packages Used**
Pathway scores were calculated in Python 3.9.17 with the following packages:
- pandas
- scipy
Plotting was performed in R 4.2.1 using the following packages:
- dplyr
- ggplot2
- corrplot

## Supplementary Figure: Odds Ratio Plot
**Overview**
Five pre-assembled lists of pathways are assembled into a dictionary mapping pathway to a list of genes in the pathway. Previously-assembled lists of likely loss-of-function (lof), gain-of-function (gof), germline, and somatic variants are imported. For each variant list (gof, lof, somatic, germline), a hypergeometric distribution is run measuring the intersection of genes in a pathway and the variant list. All p-values are FDR corrected and assembled into a single dataframe. Fisher's exact tests are run for each combination of variant list, quantifying the enrichment of gof/lof variants in the lists of germline/somatic variants. Results from the Fisher's exact tests are plotted as heatmaps.

**Packages Used**
The following packages were used in Python 3.9.17:
- pandas
- scipy
- statsmodels
- matplotlib
- seaborn
