import pandas as pd
from scipy.stats import zscore

# reformat the pathways dataframe, turn into dictionary of pathway:gene list
pathways = pd.read_csv("/wistar/auslander/McKenna/pdx/data/mel_pathway/selected_pathways.csv")
pathways = pathways.set_index('pathway_name')
pathways.index.name = None
path_dict = {}
all_genes = set()

for index, row in pathways.iterrows():
    non_na_values = [value for value in row if pd.notna(value)]
    path_dict[index] = non_na_values
    all_genes.update(set(non_na_values))


############################################################
#                  whole exome sequencing                  #  
############################################################
wes = pd.read_csv("/wistar/auslander/melanoma_pdx_m/omics_project/pdx_data_for_paper/wes.csv")
wes_copy = wes.copy() # just in case
wes.set_index("PDX", inplace=True)
wes.index.name = None
genes = [col.split('_')[0] for col in wes.columns]
wes.columns = genes

# combine columns of the same gene
wg = pd.DataFrame()
for gene in set(genes):
    wg[gene] = wes.filter(regex=f'^{gene}').max(axis=1)

wg.to_csv("/wistar/auslander/McKenna/pdx/data/mel_pathway/for_paper/wes_cols_combined.csv")

wes_cols = set(wg.columns).intersection(all_genes)
wg = wg[wes_cols]

# for each pathway, get the genes and count the number of nonzero values
nm = pd.DataFrame(index=wg.index, columns=path_dict.keys())
for id in wg.index:
    # loop through the pathway dictionary
    for pathway, genes in path_dict.items():
        valid_genes = [gene for gene in genes if gene in wg.columns]
        nm.loc[id, pathway] = wg.loc[id, valid_genes].sum()

nm.to_csv("/wistar/auslander/McKenna/pdx/outputs/mel_pathway/for_paper/wes_num_mut.csv")
nm_numeric = nm.apply(pd.to_numeric, errors='coerce')
nm_z = pd.DataFrame(nm_numeric.apply(zscore))
nm_z.to_csv("/wistar/auslander/McKenna/pdx/outputs/mel_pathway/for_paper/wes_zscored_num_mut.csv")

############################################################
#                        rnaseq                            #  
############################################################
rnaseq_raw = pd.read_csv("/wistar/auslander/melanoma_pdx_m/omics_project/pdx_data_for_paper/rnaseq.csv")
rnaseq_raw.rename(columns={rnaseq_raw.columns[0]: "PDX"}, inplace = True)
rnaseq_raw.set_index('PDX', inplace=True)
rnaseq = pd.DataFrame(rnaseq_raw.apply(zscore))
rnaseq.index.name = None

zscore_means = pd.DataFrame(index=rnaseq.index, columns=path_dict.keys())
for id in rnaseq.index:
    for pathway, genes in path_dict.items():
        valid_genes = [gene for gene in genes if gene in rnaseq.columns]
        zscore_means.loc[id, pathway] = rnaseq.loc[id, valid_genes].mean()    
        
zm_numeric = zscore_means.apply(pd.to_numeric, errors='coerce')
zm_numeric.to_csv("/wistar/auslander/McKenna/pdx/outputs/mel_pathway/for_paper/rnaseq_zscore_means.csv")

############################################################
#                          rppa                            #  
############################################################
rppa_norm_gene = pd.read_csv("/wistar/auslander/melanoma_pdx_m/omics_project/pdx_data_for_paper/rppa_norm_gene.csv")
rppa_norm_gene.set_index("Line", inplace=True)
rppa_norm_gene.index.name = None
rppa_norm_gene = rppa_norm_gene.drop(columns='Batch')
rppa_norm_gene = pd.DataFrame(rppa_norm_gene.apply(zscore))

rppa_means = pd.DataFrame(index=rppa_norm_gene.index, columns=path_dict.keys())
for id in rppa_means.index:
    for pathway, genes in path_dict.items():
        # only include the genes from a pathway that have RPPA data
        valid_genes = [gene for gene in genes if gene in rppa_norm_gene.columns]
        rppa_means.loc[id, pathway] = rppa_norm_gene.loc[id, valid_genes].mean()
        

rppa_m_numeric = rppa_means.apply(pd.to_numeric, errors='coerce')
rppa_m_numeric.to_csv("/wistar/auslander/McKenna/pdx/outputs/mel_pathway/for_paper/rppa_gene_means.csv")



############################################################
#                rppa (averaging across isoforms)          #  
############################################################
rppa_norm_gene = pd.read_csv("/wistar/auslander/melanoma_pdx_m/omics_project/pdx_data_for_paper/rppa_norm_gene.csv")
rppa_norm_gene.set_index("Line", inplace=True)
rppa_norm_gene.index.name = None
rppa_norm_gene = rppa_norm_gene.drop(columns='Batch')
rppa_norm_gene = rppa_norm_gene.groupby(rppa_norm_gene.columns.str.split('.').str[0], axis=1).mean()
rppa_norm_gene = pd.DataFrame(rppa_norm_gene.apply(zscore))

rppa_means = pd.DataFrame(index=rppa_norm_gene.index, columns=path_dict.keys())
for id in rppa_means.index:
    for pathway, genes in path_dict.items():
        # only include the genes from a pathway that have RPPA data
        valid_genes = [gene for gene in genes if gene in rppa_norm_gene.columns]
        rppa_means.loc[id, pathway] = rppa_norm_gene.loc[id, valid_genes].mean()
        

rppa_m_numeric = rppa_means.apply(pd.to_numeric, errors='coerce')
rppa_m_numeric.to_csv("/wistar/auslander/McKenna/pdx/outputs/mel_pathway/for_paper/rppa_isogene_means.csv")