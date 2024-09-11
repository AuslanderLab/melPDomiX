import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import pandas as pd
import seaborn as sns
import gseapy as gp


plt.rcParams["font.family"] = "Arial"
mpl.rcParams['pdf.fonttype'] = 42

##########################################################
#Pathway analysis for each gene on aggregated variant score

selected_pathways = pd.read_csv('selected_pathways.csv', header=None)
pathway_dict = {}
for _, row in selected_pathways.iterrows():
    pathway = row[0]
    genes = row[1:].dropna().tolist()
    pathway_dict[pathway] = genes
dd = pd.read_csv('variant_prediction_by_fa_rppa.csv')

gof_lof_genes = dd['variants'].apply(lambda x: x.split('_')[0]).unique().tolist()
selected_genes_ora_temp = [gene for genes in pathway_dict.values() for gene in genes]
ora_results_temp = gp.enrich(gene_list=gof_lof_genes,
                             gene_sets=pathway_dict,
                             background=list(selected_genes_ora_temp),
                             outdir=None)

#Filtered pathways for those with more than 1 gene in overlap
pathway_dict_filt = {key: value for key, value in pathway_dict.items() if key not in (list(ora_results_temp.results.Term.iloc[[6,8,22]])+["GO_METALLOPEPTIDASE_ACTIVITY"])}
selected_genes_ora = [gene for genes in pathway_dict_filt.values() for gene in genes]
background_genes_ora = set(selected_genes_ora).union(set(gof_lof_genes))
ora_results = gp.enrich(gene_list=gof_lof_genes,
                             gene_sets=pathway_dict_filt,
                             background=list(background_genes_ora),
                             outdir=None)

def extract_ora_pvalues(ora_results):
    pvals = pd.DataFrame(index=ora_results.res2d['Term'])
    pvals['pval'] = list(ora_results.res2d['Adjusted P-value'])
    return pvals

ora_pvals = extract_ora_pvalues(ora_results)
ora_results.results['-log10(pval)'] = -np.log10(ora_results.results['Adjusted P-value'])
ora_results.results['% Genes in set'] = ora_results.results['Overlap'].apply(lambda x: int(x.split('/')[0]) / int(x.split('/')[1]))

ora_results_plot = ora_results.results[['Term', '-log10(pval)', '% Genes in set', 'Odds Ratio']]
ora_results_plot = ora_results_plot.sort_values(by='-log10(pval)', ascending=False)
ora_results_plot['Size Category'] = pd.cut(ora_results_plot['% Genes in set'],
                                           bins=[0, 0.05, 1],
                                           labels=['< 5%', '>= 5%'])

##########################################################
#Plot significant pathways
plt.figure(figsize=(12, 8))
norm = plt.Normalize(vmin=ora_results_plot['Odds Ratio'].min(), vmax=ora_results_plot['Odds Ratio'].max())
sm = plt.cm.ScalarMappable(cmap=sns.color_palette('blend:blue,red', as_cmap=True) , norm=norm)
ax = sns.scatterplot(data=ora_results_plot, x='-log10(pval)', y='Term', hue='Odds Ratio', palette=sns.color_palette('blend:blue,red', as_cmap=True), size='Size Category', sizes={'< 5%': 50, '>= 5%': 120}, legend='full')
ax.axvline(x=-np.log10(0.05), color='black', linestyle='--')
ax.set_title('Pathways Enriched in Driver Genes with GOF/LOF Variants')
ax.set_xlabel('-log10(p-value)')
ax.set_ylabel('Pathway')
ax.get_legend().remove()
cbar = ax.figure.colorbar(sm, ax=ax)
cbar.set_label('Odds Ratio', rotation=270)
ax.legend(ax.get_legend_handles_labels()[0][-3:], ax.get_legend_handles_labels()[1][-3:], title='% Genes in set', bbox_to_anchor=(1.2, 1),borderaxespad=0.)
plt.tight_layout()
plt.show()

#Harmonize variants to genes
variant_scores = dd
pathway_gene_scores = {}
for pathway, genes in pathway_dict_filt.items():
    pathway_gene_scores[pathway] = {}
    for gene in genes:
        variant_scores_for_gene = variant_scores[variant_scores['variants'].str.startswith(gene + '_')]
        if not variant_scores_for_gene.empty:
            avg_score = variant_scores_for_gene['score'].mean()
            pathway_gene_scores[pathway][gene] = avg_score
        else:
            pathway_gene_scores[pathway][gene] = 0

#Plot heatmap for pathway analysis
heatmap_data = pd.DataFrame(pathway_gene_scores).T.fillna(0)
heatmap_data = heatmap_data[list(set(gof_lof_genes)&set(selected_genes_ora))]
heatmap_data.drop(["GO_METALLOPEPTIDASE_ACTIVITY"], axis=0, inplace=True) #no overlap
heatmap_data = heatmap_data.loc[ora_results_plot.Term]

col_linkage = sns.clustermap(heatmap_data, center=0, vmin=-0.3, vmax=0.3,cmap="coolwarm", linewidths=.5, col_cluster=True, row_cluster=False)
col_order = col_linkage.dendrogram_col.reordered_ind
heatmap_data = heatmap_data.iloc[:, col_order]

plt.figure(figsize=(15, 10))
g1p = sns.heatmap(heatmap_data, vmin=-0.3, vmax=0.3, center=0, cmap="coolwarm", cbar_kws={'label': 'Average Score'}, linewidths=.5)
plt.title('Average Score of Variants in Pathways')
plt.xlabel('Genes')
plt.ylabel('Pathways')
plt.xticks(rotation=90)
plt.show()

##########################################################
#Adjacency matrix for network visualization in Cytoscape
variants = dd['variants'].unique()
adj_matrix = pd.DataFrame(0, index=variants, columns=variants)
final_class_dict = dd.set_index('variants')['class'].to_dict()
factors = ff['Factor'].unique()
for factor in factors:
    factor_data = ff[ff['Factor'] == factor]
    for i, var1 in enumerate(variants):
        for j, var2 in enumerate(variants):
            if var1 != var2 and var1 in factor_data['Variant'].values and var2 in factor_data['Variant'].values:
                var1_label = factor_data[factor_data['Variant'] == var1]['Label'].values[0]
                var2_label = factor_data[factor_data['Variant'] == var2]['Label'].values[0]
                if (final_class_dict[var1] == var1_label) and (final_class_dict[var2] == var2_label):
                    adj_matrix.at[var1, var2] += 1
adj_matrix.to_csv("adj_matrix_factors.csv")