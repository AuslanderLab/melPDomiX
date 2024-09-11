import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import roc_auc_score
from sklearn.metrics import roc_curve, auc, precision_recall_curve, average_precision_score


plt.rcParams["font.family"] = "Arial"
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['axes.facecolor'] = 'none'

##########################################################
#Variant classification across omics using MOFA factor loadings

#MOFA output from 01_mofa_loadings
ff = pd.read_csv('gof_lof_variants_weights_scaled_all.csv')
ff.index = ff.Variant

#True GOF/LOF status for experimentally validated variants
vars = pd.DataFrame({"variants":
                     ['ISX_p.R86C','APOB_p.R1388C','TP53_p.R213*','MACF1_p.C847F','DNAJC5B_p.E22K','CHST15_p.E66K','ALK_p.K1491R','THSD7B_p.P1116S','KEAP1_p.C115R','ATXN3_p.G306R','TRHDE_p.R928Q','NIPBL_p.N674S','SYT1_p.P180S','DNAH3_p.E3387K','FCAR_p.S269G','ALK_p.D1529E','PIK3R2_p.A727T','PTCHD3_p.S309fs','KMT2B_p.P1829L','CPXM2_p.E538K','TEX15_p.E511K','TM9SF1_p.R215H','NRAS_p.Q61R','CDKN2A_p.R58*','FREM2_p.R1166C','USH2A_p.E4671K','CACNA2D3_p.E912K','DPPA3_p.S23F','BRCA1_p.D693N','KDM6A_p.T754K','RAC1_p.P29S','NRAS_p.Q61K','SLC9A4_p.G580E','MSH3_p.Q949R','SCN2A_p.R379C','OR6C1_p.S70L','SLC9C2_p.G200E','ATM_p.D1853N','BRAF_p.V600E','ADAM18_p.S244F','BRAF_p.V600K','CDKN2A_p.P114L','BRCA1_p.Q356R','PCDHAC1_p.E451K'],
                     "class":['LOF','LOF','LOF','LOF','LOF','LOF','LOF','LOF','LOF','LOF','LOF','LOF','LOF','LOF','LOF','LOF','LOF','LOF','LOF','LOF','LOF','GOF','GOF','GOF','GOF','GOF','GOF','GOF','GOF','GOF','GOF','GOF','GOF','GOF','GOF','GOF','GOF','GOF','GOF','GOF','GOF','GOF','GOF','GOF']})


sgof = set(vars['variants'][vars['class']=='GOF'])&set(ff['Variant'])
slof = set(vars['variants'][vars['class']=='LOF'])&set(ff['Variant'])

rppa_weights = [str(x).split(";") for x in ff["RPPA_Weights"]]
ff["RPPA_mean"] = [np.mean([float(l) for l in sublist] or 0) for sublist in rppa_weights]
ff["RPPA_mean"] = ff["RPPA_mean"].fillna(0)

combined_score = ff["Variant_Weight"]*ff["Gene_Weight"] + ff["Variant_Weight"]*ff["RPPA_mean"] + ff["Gene_Weight"]*ff["RPPA_mean"]
ff['mlt'] = combined_score/max(abs(combined_score))   
    
uvars = list(set(ff.index))

var_class = []
var_score = []
for v in range(len(uvars)):

    sorted_vrs = ff.loc[uvars[v]].mlt.sort_values()

    if sorted_vrs.min() < 0 and sorted_vrs.max()>0: #if variant classified as both LOF and GOF initially by MOFA
        if abs(sorted_vrs.iloc[0])>abs(sorted_vrs.iloc[-1]):
            var_class.append('LOF')
            var_score.append(sorted_vrs.iloc[0])
        else:
            var_class.append('GOF')
            var_score.append(sorted_vrs.iloc[-1])
    if sorted_vrs.min() < 0  and sorted_vrs.max()<0: #variant classified only as LOF
        var_class.append('LOF')
        var_score.append(sorted_vrs.iloc[0])
    if sorted_vrs.min() >= 0  and sorted_vrs.max()>=0: #variant classified only as GOF
        var_class.append('GOF')
        var_score.append(sorted_vrs.iloc[-1])


dd = pd.DataFrame({'variants':uvars,'class':var_class,'score':var_score})
dd.index = dd['variants']
dd["Factor"] = [ff.loc[var][ff.loc[var]["mlt"] == dd.loc[var]["score"]]["Factor"][0] for var in dd.index]
dd.to_csv('variant_prediction_by_fa_rppa.csv')


##########################################################
#Plot scores
top_gof = dd[dd['class'] == 'GOF'].nlargest(8, 'score')
top_lof = dd[dd['class'] == 'LOF'].nsmallest(8, 'score')
top_variants = pd.concat([top_gof, top_lof])
top_variants['abs_score'] = top_variants['score'].abs()
top_variants_gof = top_variants[top_variants['class'] == 'GOF'].sort_values(by='score', ascending=False)
top_variants_lof = top_variants[top_variants['class'] == 'LOF'].sort_values(by='score', ascending=True)
final_variants = pd.concat([top_variants_gof, top_variants_lof])

plt.figure(figsize=(10, 8))
sns.barplot(x='score', y='variants', data=final_variants, 
            palette=['red' if x == 'GOF' else 'blue' for x in final_variants['class']])
plt.axvline(0, color='black', linewidth=0.8)  # Vertical line at 0
plt.title('Top GOF and LOF Variants by Score')
plt.xlabel('Score')
plt.ylabel('Variants')
max_score = max(final_variants['score'].abs())
plt.xlim(-max_score, max_score)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.tight_layout()
plt.show()

#Plot heatmap
top_gof = dd[dd['class'] == 'GOF'].nlargest(8, 'score')
top_lof = dd[dd['class'] == 'LOF'].nsmallest(8, 'score')
top_variants = pd.concat([top_gof, top_lof])
heatmap_data = pd.DataFrame(index=top_variants.index, columns=['WES', 'RNAseq', 'RPPA'])
for variant in top_variants.index:
    factor = top_variants.loc[variant, 'Factor']
    factor_l = ff[ff.Factor==factor]
    if len(factor_l.loc[variant]) == 2:
        rna_weight = factor_l.loc[variant].iloc[0]['Gene_Weight']
        wes_weight = factor_l.loc[variant].iloc[0]['Variant_Weight']
        rppa_weights = factor_l.loc[variant].iloc[0]['RPPA_Weights']
    else:
        rna_weight = factor_l.loc[variant]['Gene_Weight']
        wes_weight = factor_l.loc[variant]['Variant_Weight']
        rppa_weights = factor_l.loc[variant]['RPPA_Weights']
        
    if isinstance(rppa_weights, str) and ';' in rppa_weights:
        rppa_weight = np.mean([float(w) for w in rppa_weights.split(';')])
    elif pd.isnull(rppa_weights):
        rppa_weight=0
    else:
        rppa_weight = float(rppa_weights)
    
    heatmap_data.loc[variant] = [wes_weight, rna_weight, rppa_weight]
heatmap_data = heatmap_data.astype(float)

plt.figure(figsize=(10, 8))
sns.heatmap(heatmap_data, cmap=sns.color_palette('coolwarm', as_cmap=True), center=0, linewidths=.5, cbar_kws={'label': 'Factor Loading'})
plt.title('Factor Loadings for Top GOF/LOF Variants Across Omics')
plt.xlabel('Omics Modality')
plt.ylabel('Variants')
plt.tight_layout()
plt.show()

##############################################
#AUROC
gof_prediction = dd.loc[set(sgof)&set(dd.index)]
lof_prediction = dd.loc[set(slof)&set(dd.index)]

labels=[1 for i in range(len(gof_prediction))]+[0 for i in range(len(lof_prediction))]
scores = list(gof_prediction['score'])+list(lof_prediction['score'])

fpr, tpr, _ = roc_curve(labels, scores)
roc_auc = auc(fpr, tpr)

precision, recall, _ = precision_recall_curve(labels, scores)
average_precision = average_precision_score(labels, scores)
baseline = sum(labels)/len(labels)

#Plot ROC
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 8))
ax1.plot(fpr, tpr, color='darkorange', lw=2, label='ROC curve (AUC = %0.2f)' % roc_auc)
ax1.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')
ax1.set_xlim([0.0, 1.0])
ax1.set_ylim([0.0, 1.05])
ax1.set_xlabel('False Positive Rate', fontsize=14)
ax1.set_ylabel('True Positive Rate', fontsize=14)
ax1.set_title('ROC Curve', fontsize=16)
ax1.legend(loc="lower right", fontsize=12)
ax1.tick_params(axis='both', which='major', labelsize=12)

ax2.plot(recall, precision, color='b', lw=2, label='Precision-Recall curve (AUPRC = %0.2f)' % average_precision)
ax2.plot([0, 1], [baseline, baseline], color='navy', lw=2, linestyle='--', label='Baseline (P/N = %0.2f)' % baseline)
ax2.set_xlim([0.0, 1.0])
ax2.set_ylim([0.0, 1.05])
ax2.set_xlabel('Recall', fontsize=14)
ax2.set_ylabel('Precision', fontsize=14)
ax2.set_title('Precision-Recall Curve', fontsize=16)
ax2.legend(loc="lower left", fontsize=12)
ax2.tick_params(axis='both', which='major', labelsize=12)

plt.tight_layout()
plt.show()