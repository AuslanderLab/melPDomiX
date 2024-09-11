import pandas as pd
import numpy as np
import pickle
from scipy.stats import hypergeom
from statsmodels.stats.multitest import fdrcorrection
from scipy.stats import fisher_exact
import matplotlib.pyplot as plt
from statsmodels.graphics.mosaicplot import mosaic
import seaborn as sns

####################################################
#            load in pathways                      #  
####################################################
# Pathway dfs
canonical = pd.read_csv("/wistar/auslander/scripts/pathways/canonical_pathways_gene_sets.txt", 
                        sep='\t', names=range(1465), low_memory=False)
canonical = canonical.drop(canonical.columns[1], axis=1)
curated = pd.read_csv("/wistar/auslander/scripts/pathways/curated_gene_sets.txt", sep = '\t', 
                      names=range(1961), header = None)
curated = curated.drop(curated.columns[1], axis = 1)
go = pd.read_csv("/wistar/auslander/scripts/pathways/GO.csv", low_memory= False)
go = go.drop(columns=['link'])
hallmark = pd.read_csv("/wistar/auslander/scripts/pathways/hallmark.txt", 
                       sep='\t', header = None)
hallmark = hallmark.drop(hallmark.columns[1], axis=1)
onco = pd.read_csv("/wistar/auslander/scripts/pathways/onco_pathways.csv")
onco = onco.transpose()
onco.reset_index(inplace=True)

####################################################
#     load in lists from supplemental dataframe    #  
####################################################
variants = pd.read_csv("/wistar/auslander/melanoma_pdx_m/omics_project/manuscript/supplementary/variant_classification_full.csv", index_col = 0)

gof = set(variants[variants['GOF_LOF'] == 'GOF']['variants'].tolist())
genes_gof = list(set([i.split('_')[0] for i in gof]))

lof = set(variants[variants['GOF_LOF'] == 'LOF']['variants'].tolist())
genes_lof = list(set([i.split('_')[0] for i in lof]))

germline = set(variants[variants['tcag_classification'] == 'germline']['variants'].tolist())
genes_germline = list(set([i.split('_')[0] for i in germline]))

somatic = set(variants[variants['tcag_classification'] == 'somatic']['variants'].tolist())
genes_somatic = list(set([i.split('_')[0] for i in somatic]))

####################################################
# make dictionary of all pathways and their genes  #  
####################################################
pathway_dict = {}
gene_set = set()
def make_dictionary(df):
     # iterate through each row
     for index, row in df.iterrows():
         values = []
         pathway_name = row.iloc[0]
         values = row.iloc[1:].tolist() # get a list of the genes in that pathway
         values = [v for v in values if pd.notna(v)] # remove NA values
         for value in values:
             gene_set.add(value) # add to list that keeps track of total num genes
         pathway_dict[pathway_name] = values

# run this on all 
make_dictionary(canonical)
make_dictionary(curated)
make_dictionary(go)
make_dictionary(hallmark)
make_dictionary(onco)

####################################################
#          run hypergeometric distribution         #  
####################################################
# for each pathway:
# 20,000 – number of genes
# 20 – number of genes in a pathway
# 15 – number of GOF genes (can be LOF/germline/somatic variants)
# 2 – the number of gof genes in the pathway (intersect)
# 1-hypergeom.cdf(2, 20000, 20, 15)

gene_count = len(gene_set) # all genes across all pathways
def get_p_values(sublist, pathway_dict, gene_count):
     pvals = {} # dictionary with pathway: pval (only calculated if intersect >= 2)
     int_genes = {} # dictionary with pathway: genes in sublist + pathway
     sublist_length = len(sublist) # number of gof/lof/germline/somatic genes
     for pathway in pathway_dict.keys():
        pathway_list = pathway_dict[pathway] # get genes in that pathway
        pathway_count = len(pathway_list) # get number of genes in the pathway
        intersection_genes = list(set(sublist) & set(pathway_list)) # genes in the sublist and the pathway
        int_genes[pathway] = intersection_genes # save as a key:value pair to dictionary
        pathway_intersect = len(intersection_genes) # genes in (pathway + gof/lof/somatic/germline)
        if pathway_intersect >= 2:
            pval = 1 - hypergeom.cdf(pathway_intersect - 1, gene_count, pathway_count, sublist_length)
            pvals[pathway] = pval
     new_df = pd.DataFrame(list(pvals.items()), columns=['pathway', 'pval'])
     return new_df, int_genes


gof_p_vals, gof_int_genes = get_p_values(genes_gof, pathway_dict, gene_count)
gof_p_vals['origin'] = "gof"
gof_p_vals['genes'] = gof_p_vals['pathway'].map(gof_int_genes)

lof_p_vals, lof_int_genes = get_p_values(genes_lof, pathway_dict, gene_count)
lof_p_vals['origin'] = "lof"
lof_p_vals['genes'] = lof_p_vals['pathway'].map(lof_int_genes)

somatic_p_vals, somatic_int_genes = get_p_values(genes_somatic, pathway_dict, gene_count)
somatic_p_vals['origin'] = "somatic"
somatic_p_vals['genes'] = somatic_p_vals['pathway'].map(somatic_int_genes)

germline_p_vals, germline_int_genes = get_p_values(genes_germline, pathway_dict, gene_count)
germline_p_vals['origin'] = "germline"
germline_p_vals['genes'] = germline_p_vals['pathway'].map(germline_int_genes)

dfs = [gof_p_vals, lof_p_vals, somatic_p_vals, germline_p_vals]
all_pvals = pd.concat(dfs)

# Correct the p values --------------------------------------------------------------------------------------
def correct_p_vals(pval_df, alpha):
    pvals_list = pval_df['pval']
    rejected, corrected_p_values = fdrcorrection(pvals_list, alpha = alpha)
    pval_df['pval_corrected'] = corrected_p_values
    pval_df = pval_df.sort_values(by='pval_corrected', ascending = True)
    pval_df = pval_df.reindex(columns=['pathway', 'pval', 'pval_corrected', 'origin', 'genes'])
    # pval_df.to_csv(savepath)
    return pval_df, corrected_p_values

all_pvals_corrected, pvals = correct_p_vals(all_pvals, 0.05)

########################################################
# make contingency table based on sublist abundance    #
########################################################
def get_sublist(df, desired, undesired):
    sublist = set()
    # only pathways in desired list and not in undesired list
    pathways_desired = set(df[df['origin'] == desired]['pathway'])
    pathways_undesired = set(df[df['origin'] == undesired]['pathway'])
    pathways_only_desired = set(pathways_desired) - set(pathways_undesired)
    sublist.update(pathways_only_desired)
    return sublist

germline_list = get_sublist(all_pvals_corrected, "germline", "somatic")
somatic_list = get_sublist(all_pvals_corrected, "somatic", "germline")
gof_list = get_sublist(all_pvals_corrected, "gof", "lof")
lof_list = get_sublist(all_pvals_corrected, "lof", "gof")

# running fisher tests ---------------------------------------------------------
# contingency table
#           GOF        LOF 
#somatic     a          b
#germline    c          d

somatic_gof = len(gof_list & somatic_list)
somatic_lof = len(lof_list & somatic_list)
germline_gof = len(gof_list & germline_list)
germline_lof = len(lof_list & germline_list)

odds_ratio_pathway, pval_pathway = fisher_exact([[somatic_gof, somatic_lof], 
                                [germline_gof, germline_lof]])
contingency_table_pathway = {
    'germline/somatic': ['germline', 'somatic'],
    'GOF': [somatic_gof, somatic_lof],
    'LOF': [germline_gof, germline_lof],
}
contingency_df_pathway = pd.DataFrame(contingency_table_pathway)
contingency_df_pathway.set_index('germline/somatic', inplace=True)
contingency_df_pathway.to_csv("/wistar/auslander/McKenna/pdx/outputs/mel_pathway/for_paper/pathway_contingency.csv")

# on gene level
# contingency table
#           GOF        LOF 
#somatic     a          b
#germline    c          d

vars_somatic_gof = len(gof & somatic)
vars_somatic_lof = len(lof & somatic)
vars_germline_gof = len(gof & germline)
vars_germline_lof = len(lof & germline)

odds_ratio_var, pval_var = fisher_exact([[vars_somatic_gof, vars_somatic_lof], 
                                [vars_germline_gof, vars_germline_lof]])

contingency_table_variant = {
    'germline/somatic': ['somatic', 'germline'],
    'GOF': [vars_somatic_gof, vars_germline_gof],
    'LOF': [vars_somatic_lof, vars_germline_lof],
}
contingency_df_variant = pd.DataFrame(contingency_table_variant)
contingency_df_variant.set_index('germline/somatic', inplace=True)
contingency_df_variant.to_csv("/wistar/auslander/McKenna/pdx/outputs/mel_pathway/for_paper/variant_contingency.csv")

########################################################
#     visualize the contingency tables graphically     #
########################################################
# Pathways ---------------------------------------------------------------------
plt.clf()
sns.heatmap(contingency_df_pathway, annot=True, fmt="d", cmap="YlGnBu")
plt.title("Overlap of Pathways in Each Sublist", y = 1)
subtitle = "Odds Ratio = " + str(odds_ratio_pathway) + ", p = " + str(pval_pathway)
plt.ylabel("")
plt.xticks(ticks=[0.5, 1.5], labels=['Gain of Function', 'Loss of Function'], rotation=0)
plt.xlabel(subtitle)
plt.yticks(ticks=[0.5, 1.5], labels=['Somatic', 'Germline'], rotation=90)
plt.tight_layout(rect=[0, 0, 1, 0.90]) 
plt.savefig("/wistar/auslander/McKenna/pdx/outputs/mel_pathway/for_paper/contingency_heatmap_pathway.pdf")

# Variants ---------------------------------------------------------------------
plt.clf()
sns.heatmap(contingency_df_variant, annot=True, fmt="d", cmap="YlGnBu")
plt.title("Overlap of Variants in Each Sublist", y = 1)
subtitle = "Odds Ratio = " + str(odds_ratio_var) + ", p = " + str(pval_var)
plt.ylabel("")
plt.xticks(ticks=[0.5, 1.5], labels=['Gain of Function', 'Loss of Function'], rotation=0)
plt.xlabel(subtitle)
plt.yticks(ticks=[0.5, 1.5], labels=['Somatic', 'Germline'], rotation=90)
plt.tight_layout(rect=[0, 0, 1, 0.90]) 
plt.savefig("/wistar/auslander/McKenna/pdx/outputs/mel_pathway/for_paper/contingency_heatmap_var.pdf")