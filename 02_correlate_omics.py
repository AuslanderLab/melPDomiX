import pandas as pd
import numpy as np

############################################################
#                      load in data                        #
############################################################
wes = pd.read_csv("/wistar/auslander/McKenna/pdx/outputs/mel_pathway/for_paper/wes_zscored_num_mut.csv")
wes = wes.set_index("Unnamed: 0")
wes.index.name = None

rnaseq = pd.read_csv("/wistar/auslander/McKenna/pdx/outputs/mel_pathway/for_paper/rnaseq_zscore_means.csv")
rnaseq = rnaseq.set_index("Unnamed: 0")
rnaseq.index.name = None

rppa = pd.read_csv("/wistar/auslander/McKenna/pdx/outputs/mel_pathway/for_paper/rppa_gene_means.csv")
rppa = rppa.set_index("Unnamed: 0")
rppa.index.name = None

isorppa = pd.read_csv("/wistar/auslander/McKenna/pdx/outputs/mel_pathway/for_paper/rppa_isogene_means.csv", index_col = 0)
isorppa.index.name = None

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


# add info about omics type to each df
wes = wes.rename(columns=lambda x: x + '_wes')
rnaseq = rnaseq.rename(columns=lambda x: x + '_rnaseq')
rppa = rppa.rename(columns=lambda x: x + '_rppa')
isorppa = isorppa.rename(columns=lambda x: x + '_rppa')

omics_all = pd.concat([wes, rnaseq, rppa], axis=1)
omics_corr = omics_all.corr(method='spearman')

############################################################
#                  incorporate clinical data               # 
#                  PR = partial response                   #
#                  CR = complete response                  #
#                  PD = progressive disease                #
#                  SD = stable disease                     #
############################################################

# Load in data and reformat
clin = pd.read_csv("/wistar/auslander/melanoma_pdx_m/omics_project/pdx_data_for_paper/clin.csv")
clin = clin.set_index("WM_number")
clin.index.name = None

def create_corr_df(wes, rnaseq, rppa, treatment_type, clinical_df, resp_corr_path, nonresp_corr_path, mrsd_corr_path):
    # Reformatting the data
    if treatment_type == "targeted":
        subset_treatment = pd.DataFrame(clin['TT best response pre/post'])
        subset_treatment.rename(columns={'TT best response pre/post': 'response'}, inplace=True)
        subset_treatment = subset_treatment.dropna()
    if treatment_type == "immunotherapy":
        subset_treatment = pd.DataFrame(clinical_df['IT best response pre/post'])
        subset_treatment.rename(columns={'IT best response pre/post': 'response'}, inplace=True)
        subset_treatment = subset_treatment.dropna() 
    
    # Getting patient IDs for responders, nonresponders, and MR/SD
    responders = subset_treatment[subset_treatment['response'] == 'PR/CR'].index 
    nonresponders = subset_treatment[subset_treatment['response'] == 'PD'].index
    mr_sd = subset_treatment[subset_treatment['response'].isin(['MR', 'SD'])].index
    
    wes_responders = wes.loc[responders]
    wes_nonresponders = wes.loc[nonresponders]
    wes_mrsd = wes.loc[mr_sd]
    
    rnaseq_responders = rnaseq.loc[responders]
    rnaseq_nonresponders = rnaseq.loc[nonresponders]
    rnaseq_mrsd = rnaseq.loc[mr_sd]
    
    rppa_responders = rppa.loc[responders]
    rppa_nonresponders = rppa.loc[nonresponders]
    rppa_mrsd = rppa.loc[mr_sd]
    
    # Concatenating this all together
    omics_responders = pd.concat([wes_responders, rnaseq_responders, rppa_responders], axis=1)
    omics_nonresponders = pd.concat([wes_nonresponders, rnaseq_nonresponders, rppa_nonresponders], axis=1)
    omics_mrsd = pd.concat([wes_mrsd, rnaseq_mrsd, rppa_mrsd], axis=1)
    
    # Calculating correlation matrix
    responders_corr = omics_responders.corr(method = 'spearman')
    nonresponders_corr = omics_nonresponders.corr(method='spearman')
    mrsd_corr = omics_mrsd.corr(method = 'spearman')
    
    # Save correlation matrices
    responders_corr.to_csv(resp_corr_path)
    nonresponders_corr.to_csv(nonresp_corr_path)
    mrsd_corr.to_csv(mrsd_corr_path)
    
    # Return all three correlation matrices
    return responders_corr, nonresponders_corr, mrsd_corr

responders_corr_tt, nonresponders_corr_tt, mrsd_corr_tt = create_corr_df(wes, rnaseq, rppa, 
                                                                         "targeted", clin, 
                                                                         "/wistar/auslander/McKenna/pdx/outputs/mel_pathway/for_paper/responders_corr.csv",
                                                                         "/wistar/auslander/McKenna/pdx/outputs/mel_pathway/for_paper/nonresponders_corr.csv",
                                                                         "/wistar/auslander/McKenna/pdx/outputs/mel_pathway/for_paper/mrsd_corr.csv")
responders_corr, nonresponders_corr, mrsd_corr = create_corr_df(wes, rnaseq, rppa, "immunotherapy", clin)

isoresponders_corr, isononresponders_corr, isomrsd_corr = create_corr_df(wes, rnaseq, isorppa,
                                                                         "immunotherapy", clin, 
                                                                         "/wistar/auslander/McKenna/pdx/outputs/mel_pathway/for_paper/isoresponders_corr.csv",
                                                                         "/wistar/auslander/McKenna/pdx/outputs/mel_pathway/for_paper/isononresponders_corr.csv",
                                                                         "/wistar/auslander/McKenna/pdx/outputs/mel_pathway/for_paper/isomrsd_corr.csv")

# Changing the 2D correlation matrices for responders and nonresponders into 1D matrices
def reformat_corr(df, group):
    df_1d = df.stack().reset_index()
    if group == "responder":
        df_1d.columns = ['Row', 'Column', 'corr_coeff_resp']
        group_col = 'corr_coeff_resp'
    elif group == "nonresponder":
        df_1d.columns = ['Row', 'Column', 'corr_coeff_nr']
        group_col = 'corr_coeff_nr'
    elif group == "mr":
        df_1d.columns = ['Row', 'Column', 'corr_coeff_mr']
        group_col = 'corr_coeff_mr'
    elif group == "sd":
        df_1d.columns = ['Row', 'Column', 'corr_coeff_sd']
        group_col = 'corr_coeff_sd'
    
    df_1d = df_1d[df_1d['Row'] != df_1d['Column']]
    df_1d['pair'] = df_1d.apply(lambda x: tuple(sorted([x['Row'], x['Column']])), axis=1)
    df_1d = df_1d.drop_duplicates(subset=['pair'])
    df_1d['index'] = df_1d['Row'] + ' | ' + df_1d['Column']
    df_1d.set_index('index', inplace=True)
    df_1d = df_1d[[group_col]]
    df_1d = df_1d.sort_values(by=group_col)
    df_1d.index.name = None
    
    return df_1d

def make_1d(responders_corr, nonresponders_corr, mrsd_corr, resp_nonresp_path, mrsd_path):
    responders_1d = reformat_corr(responders_corr, "responder")
    nonresponders_1d = reformat_corr(nonresponders_corr, "nonresponder")
    
    mr_1d = reformat_corr(mrsd_corr, "mr")
    sd_1d = reformat_corr(mrsd_corr, "sd")
    
    # Calculate product and difference of responder/nonresponder correlation coefficients
    resp_nonresp = pd.concat([responders_1d, nonresponders_1d], axis=1)
    resp_nonresp['product'] = resp_nonresp['corr_coeff_resp'] * resp_nonresp['corr_coeff_nr']
    resp_nonresp = resp_nonresp.sort_values(by='product', ascending=False)
    resp_nonresp['diff'] = resp_nonresp['corr_coeff_resp'] - resp_nonresp['corr_coeff_nr']
    resp_nonresp.to_csv(resp_nonresp_path)
    
    # Calculate product and difference of mixed response/stable disease correlation coefficients
    mr_sd = pd.concat([mr_1d, sd_1d], axis=1)
    mr_sd['product'] = mr_sd['corr_coeff_mr'] * mr_sd['corr_coeff_sd']
    mr_sd = mr_sd.sort_values(by='product', ascending=False)
    mr_sd['diff'] = mr_sd['corr_coeff_mr'] - mr_sd['corr_coeff_sd']
    mr_sd.to_csv(mrsd_path)

    
make_1d(responders_corr, nonresponders_corr, mrsd_corr,
        "/wistar/auslander/McKenna/pdx/outputs/mel_pathway/for_paper/resp_nonresp_pathways.csv"
        "/wistar/auslander/McKenna/pdx/outputs/mel_pathway/for_paper/mrsd_pathways.csv")