# Variant Classification

## Overview

This README describes the code used to classify variant GOF/LOF status and generate Figure 4 / Extended Figure 4.

## MOFA embedding

**01_mofa_loadings.R**

Using a concatenated dataset of batch corrected and normalized WES, RNAseq, and RPPA data (see Methods), Multi-Omics Factor Analysis (MOFA) (Argelaguet et al., 2018) is used to learn omics-specific factor loadings. Gene names are harmonized across omics, and a preliminary GOF/LOF status is derived for each variant based on concordant pairwise signs of factor loadings across omics.

**02_variant_scoring.py**

Variant weights across omics from MOFA model are used to create a score that annotates GOF/LOF status (see Methods). Descriptive visualizations for variant classifications are detailed. AUROC and precision-recall is computed for classification metric of variants with ground truth GOF/LOF label. 

**03_downstream_analysis.R**

Using a database of 31 pre-selected pathways, 26 pathways are overlapped with driver gene set used to compute MOFA embedding. Pathway enrichment for each gene is calculated by computing multiply-corrected enrichment analysis on said gene set comrpised of genes with GOF/LOF variants, and a score for each gene is computed by taking the average over all variants. Adjacency matrix for number of GOF/LOF concurrent variant statuses across factors is computed for network visualization in Cytoscape.

**Packages Used**

MOFA embedding and processing was computed in R 4.3.2 using the following packages:
- MOFA2
- dplyr
  
Scores, downstream analyses, and visualizations were computed in Python 3.11.7 with the following packages:
- pandas
- matplotlib
- seaborn
- sklearn.metrics
- gseapy
