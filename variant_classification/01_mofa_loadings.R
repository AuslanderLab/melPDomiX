library(MOFA2)
library(dplyr)

#####################################
#MOFA on combined WES, RNAseq, and RPPA data after normalization
#Take combined_zscore_data.csv() to be the concatenated omics dataset with normalized, batch corrected data
setwd()
combined_zscore_data <- read.csv("combined_zscore_data.csv")
combined_zscore_data$X <- NULL
combined_zscore_data$group <- NULL

#Using default parameters from MOFA documentation with sparse weights
MOFAobject <- create_mofa(combined_zscore_data)
print(MOFAobject)
plot_data_overview(MOFAobject)
data_opts <- get_default_data_options(MOFAobject)
head(data_opts)
model_opts <- get_default_model_options(MOFAobject)
model_opts$num_factors <- 10
model_opts$spikeslab_weights <- TRUE
model_opts$spikeslab_factors <- FALSE
head(model_opts)
model_opts$likelihoods
model_opts$likelihoods["WES"][[1]] = "bernoulli"
train_opts <- get_default_training_options(MOFAobject)
train_opts$drop_factor_threshold <- -1
head(train_opts)
MOFAobject <- prepare_mofa(
  object = MOFAobject,
  data_options = data_opts,
  model_options = model_opts,
  training_options = train_opts
)
outfile = file.path(getwd(),"model_all6.hdf5")
MOFAobject.trained <- run_mofa(MOFAobject, outfile, use_basilisk=TRUE)
model <- load_model(outfile)

#Factor loading matrices
rna_df <- get_weights(model, views = "all", factors = "all")[[1]]
rppa_df <- get_weights(model, views = "all", factors = "all")[[2]]
wes_df <- get_weights(model, views = "all", factors = "all")[[3]]

minMax <- function(x) {
  x/max(abs(x))
}
rna_df <- apply(rna_df, MARGIN = 2, function(x) minMax(x))
rppa_df <- apply(rppa_df, MARGIN = 2, function(x) minMax(x))
wes_df <- apply(wes_df, MARGIN = 2, function(x) minMax(x))
Gene <- c('YWHAB','YWHAZ','EIF4EBP1','EIF4EBP1','EIF4EBP1','TP53BP1','ARAF','ARAF','ACACA/ACACB','ACACA/ACACB','ACSL1','ACVRL1','AKT1/2/3','AKT1','AKT1','AKT2','AKT2','AKT1/2/3','AKT1/2/3','AMBRA1','PRKAA2','PRKAA1/2','PRKAA1/2','AR','ARID1A','ASNS','ATG3','ATG4B','ATG5','ATG7','ATM','ATR','ATRX','ATR','AURKA','AURKB','AXL','ACTB','CTNNB1','CTNNB1','BRAF','BRAF','CD276','VTCN1','BAD','BAK1','BAX','BCL2L1','BCL2A1','BECN1','BID','BCL2L11','MAPK7','BRD4','ABL1','ABL1','BIRC3','JUN','KIT','MET','MYC','RAF1','RAF1','CA9','CASP3','CASP7','CASP8','CAV1','MS4A1','DPP4','CD38','CD4','CD44','CD86','CDC25C','CDK1','CDC42/RAC1','CDC6','CDK1/2/3','CDKN2A','CDT1','CHEK1','CHEK1','CHEK2','CIITA','CLDN7','COG3','COL6A1','GJA1','NR2F2','COX4I1','PTGS2','CREB1','CREB1','CSK','RBBP8','CCNB1','CCND1','TUBA4A/TUBA3C','DAPK2','DDB1','DDR1','PARK7','HIST1H3A','LIG4','POLG','DNMT1','DNM1L','DUSP4','DUSP6','DVL3','CDH1','EEF2','EEF2K','EGFR','EGFR','EIF4E','EIF4E','EIF4G1','ELK1','ENO1','ENO2','EPHA2','ERS1','ESR1','ERCC5','MAPK7','ESRRA','ETS1','PTK2','PTK2','FASN','FGF2','FN1','TNFRSF12A','FOXM1','FOXO3','FOXO3','FRS2','G6PD','GAB2','GAPDH','GATA6','GCLC','GCLM','KAT2A','GLI1','GLI3','GLUD1','GLS','GZMB','GRB7','HSPA9','GSK3A/GSK3B','GSK3B','GYS1','GYS1','ERBB2','ERBB3','ERBB3','NRG1','HES1','HK1','HK2','HIST1H3A-F','ARHGAP45','HSBP1','HSPD1','HSPA1A','IDO1','IGF1R/INSR','IGFBP2','IGF1R','IL6','INPP4B','INSR','IRF1','IRS1','IRS2','JAG1','JAK2','MAPK9','MAPK8','TRIM28','LAD1','HUWE1','MAP1LC3A/B','LCK','LDHA','LRP6','MAPK1/MAPK3','MCL1','SLC16A3','MDM2','MAP2K1','MAP2K1/MAP2K1','MAP2K2','MLANA','BABAM1','BABAM1','NF2','MIF','MITF','MFN1','MFN2','MLKL','MMP14','MMP2','MKNK1','MRAP','MSH2','MSH6','MSI2','MTOR','MTOR','MYH11','MYH9','MYH9','PKMYT1','CDH2','NAPSA','NDRG1','RELA','NOTCH1','NOTCH1','NOTCH3','NFE2L2','POU5F1','CDH3','CDKN1A','CDKN1B','CDKN1B','CDKN1B','MAPK14/11/12','MAPK11/12/13/14','MAPK1/MAPK3','TP53','RPS6KB1','RPS6KB1','RPS6K','PAICS','PAK1','PAK4','PAR','PARG','PARP1','PTCH1','PXN','PDCD1','CD274','PDCD4','PDHA1','PDHK1','PDPK1','PDPK1','PEA15','PEA15','PHGDH','PHLPP1','PIK3CA','PIK3R1','PRKAR1A','PRKCA/PRKCB','PRKCA/B/D/E/H/Q','PRKCD','PRKCA','PKM','PLCG2','PLK1','PMS2','PGR','AKT1S1','PREX1','PTEN','PTPN12','BBC3','PYGB','PTK2B','RAB11A/B','RAB25','RAD23A','RAD50','RAD51','RPTOR','RBM15','RB1','RICTOR','RICTOR','RIP','RIPK3','RPA2','RPA2','RRM1','RRM2','RPS6KA1/2/3','S100A4','RPS6','RPS6','SDHA','SETD2','SFRP1','SHC1','PTPN11','PTPN11','SLC1A5','SLFN11','SMAD1','SMAD3','SOD2','SOX2','SRC','SRC','STAT3','STAT3','STAT5A','STMN1','TMEM173','WWTR1','TFAM','TFRC','TIGAR','TRIM25','TRIP13','TSC1','NKX2-1','TSC2','TSC2','TUFM','TYRO3','HIST1H2BB','UBAC1','ULK1','UVRAG','VASP','VAV1','KDR','WEE1','WEE1','WIPI1','WIPI2','XBP1','XIAP','ERCC4','XRCC1','YAP1','YAP1','YBX1','ZAP70','ACSS2','ACLY','ATM','ATP5PD','ATR','AURKA-C','BCL2','CANX','TNFRSF4','CDK9','CENPA','CHD1L','CRABP1','CRABP2','CCNE1','DDR1','DVL3','DYRK1B','E2F1','EPHA2','EPHA2','EPHA2','ERCC1','MECOM','FABP5','FANCD2','FRS2','GCLM','NR3C1','GRB2','H2AX','ERBB2','HIF1A','H3C1','HLA-DQA1','HLA-DRA','HNRNPK','INPP4B','IRF1','IRF3','MAPK9','MAPK8','KEAP1','LCN2','LYN','MACC1','MAPK1/MAPK3','MAP2K1/MAP2K2','PMEL','MDK','ERRFI1','NDUFB4','CDKN1B','MAPK11/12/13/14','PAK4-6','PAK1-3','PAX6','PAX8','EIF2AK3','PGM1','PIK3R1','PLCG1','PLCG1','PLCG2','PRC1','RAD17','RPS6KA1','SFRP1','SGK1','SGK3','SHC1','SMAD4','SOX17','SOX7','STAT1','STAT3','SYP','TEAD1-4','KDR','VHL','YES1','ZAP70','ZEB1','ALKBH5','CGAS','PTGS2','FTO','GLI3','IGF2BP3','METTL3','PIP4K2A','PIP4K2B','SRSF1','SIRPA','TRIM24','WTAP','YTHDF2','YTHDF3','FLCN','PRMT1','PRMT5','RAB11FIP1','SRC','TSC1','FGFR1','FGFR2','RORA','TACSTD2','CDH6','CASP7','EZH2','FOLR1','KDR')
RPPA <- c('14-3-3-beta-R-V','14-3-3-zeta-R-V','4E-BP1-R-V','4E-BP1_pS65-R-V','4E-BP1_pT37_T46-R-V','53BP1-R-V','A-Raf-R-V','A-Raf_pS299-R-C','ACC1-R-C','ACC_pS79-R-V','ACSL1-R-V','ACVRL1-R-C','Akt-R-V','Akt1-R-V','Akt1_pS473-R-V','Akt2-R-V','Akt2_pS474-R-C','Akt_pS473-R-V','Akt_pT308-R-V','Ambra1_pS52-R-C','AMPK-a2_pS345-R-V','AMPKa-R-C','AMPKa_pT172-R-C','AR-R-V','ARID1A-R-C','ASNS-R-V','Atg3-R-V','Atg4B-R-C','Atg5-R-C','Atg7-R-V','ATM-R-V','ATR-R-C','ATRX-R-C','ATR_pS428-R-C','Aurora-A-R-C','Aurora-B-R-V','Axl-R-V','b-Actin-R-C','b-Catenin-R-V','b-Catenin_pT41_S45-R-V','B-Raf-R-C','B-Raf_pS445-R-V','B7-H3-R-C','B7-H4-R-C','Bad_pS112-R-V','Bak-R-C','Bax-R-V','Bcl-xL-R-V','BCL2A1-R-V','Beclin-R-C','Bid-R-C','Bim-R-V','BMK1-Erk5_pT218_Y220-R-V','BRD4-R-V','c-Abl-R-V','c-Abl_pY412-R-C','c-IAP2-R-C','c-Jun_pS73-R-V','c-Kit-R-V','c-Met_pY1234_Y1235-R-V','c-Myc-R-C','C-Raf-R-C','C-Raf_pS338-R-V','CA9-R-C','Caspase-3-cleaved-R-C','Caspase-7-cleaved--R-C','Caspase-8-cleaved-R-C','Caveolin-1-R-V','CD20-R-C','CD26-R-V','CD38-R-C','CD4-R-V','CD44-R-C','CD86-R-C','cdc25C-R-V','cdc2_pY15-R-C','Cdc42-R-C','Cdc6-R-V','CDK1_pT14-R-C','CDKN2A-R-C','CDT1-R-V','Chk1_pS296-R-V','Chk1_pS345-R-C','Chk2_pT68-R-C','CIITA-R-C','Claudin-7-R-V','COG3-R-V','Collagen-VI-R-V','Connexin-43-R-C','Coup-TFII-R-C','Cox-IV-R-V','Cox2-R-C','Creb-R-C','CREB_pS133-R-C','CSK-R-C','CtIP-R-V','Cyclin-B1-R-V','Cyclin-D1-R-C','D-a-Tubulin-R-V','DAPK2-R-C','DDB-1-R-V','DDR1-R-V','DJ1-R-V','DM-Histone-H3-R-V','DNA-Ligase-IV-R-C','DNA_POLG-R-V','DNMT1-R-V','DRP1-R-V','DUSP4-R-V','DUSP6-R-C','DVL3-R-V','E-Cadherin-R-V','eEF2-R-C','eEF2K-R-V','EGFR-R-V','EGFR_pY1173-R-V','eIF4E-R-V','eIF4E_pS209-R-V','eIF4G-R-C','Elk1_pS383-R-C','Enolase-1-R-V','Enolase-2-R-V','EPHA2-R-C','ER-a-R-V','ER-a_pS118-R-V','ERCC5-R-C','Erk5-R-V','ERRalpha-R-V','Ets-1-R-V','FAK-R-C','FAK_pY397-R-V','FASN-R-V','FGF-basic-R-C','Fibronectin-R-V','FN14-R-C','FOXM1-R-V','FOXO3-R-V','FoxO3a_pS318_S321-R-C','FRS2-a_pY196-R-C','G6PD-R-V','Gab2-R-V','GAPDH-R-C','GATA6-R-V','GCLC-R-C','GCLM-R-C','GCN5L2-R-V','Gli1-R-C','Gli3-R-C','Glutamate-D1-2-R-V','Glutaminase-R-C','Granzyme-B-R-V','GRB7-R-V','Grp75-R-C','GSK-3a-b_pS21_S9-R-V','GSK-3B-R-C','Gys-R-V','Gys_pS641-R-V','HER2_pY1248-R-C','HER3-R-V','HER3_pY1289-R-C','Heregulin-R-V','HES1-R-V','Hexokinase-I-R-C','Hexokinase-II-R-V','Histone-H3-R-V','HMHA1-R-V','HSP27_pS82-R-V','HSP60-R-V','HSP70-R-C','IDO-R-C','IGF1R_pY1135_Y1136-R-V','IGFBP2-R-V','IGFRb-R-C','IL-6-R-C','INPP4b-R-V','IR-b-R-C','IRF-1-R-C','IRS1-R-V','IRS2-R-C','Jagged1-R-V','Jak2-R-V','JNK2-R-C','JNK_pT183_Y185-R-V','KAP1-R-V','LAD1-R-V','Lasu1-R-V','LC3A-B-R-C','Lck-R-V','LDHA-R-C','LRP6_pS1490-R-V','MAPK_pT202_Y204-R-V','Mcl-1-R-V','MCT4-R-V','MDM2_pS166-R-V','MEK1-R-V','MEK1_p_S217_S221-R-V','MEK2-R-V','MelanA-R-C','MERIT40-R-C','MERIT40_pS29-R-V','Merlin-R-C','MIF-R-C','MITF-R-V','Mitofusin-1-R-V','Mitofusin-2-R-V','MLKL-R-V','MMP14-R-V','MMP2-R-V','Mnk1-R-V','MRAP-R-C','MSH2-R-C','MSH6-R-C','MSI2-R-C','mTOR-R-V','mTOR_pS2448-R-C','MYH11-R-C','Myosin-IIa-R-C','Myosin-IIa_pS1943-R-V','Myt1-R-C','N-Cadherin-R-V','NAPSIN-A-R-C','NDRG1_pT346-R-V','NF-kB-p65_pS536-R-C','Notch1-R-V','Notch1-cleaved-R-V','Notch3-R-C','NRF2-R-C','Oct-4-R-C','P-Cadherin-R-C','p21-R-C','p27-Kip1-R-V','p27_pT157-R-C','p27_pT198-R-V','p38-MAPK-R-V','p38_pT180_Y182-R-V','p44-42-MAPK-R-V','p53-R-C','p70-S6K1-R-V','p70-S6K_pT389-R-V','p90RSK_pT573-R-C','PAICS-R-C','PAK1-R-V','PAK4-R-V','PAR-R-C','PARG-R-C','PARP-R-V','Patched-R-C','Paxillin-R-C','PD-1-R-V','PD-L1-R-C','Pdcd4-R-C','PDHA1-R-V','PDHK1-R-C','PDK1-R-V','PDK1_pS241-R-V','PEA-15-R-V','PEA-15_pS116-R-V','PHGDH-R-C','PHLPP-R-V','PI3K-p110-a-R-C','PI3K-p85-R-V','PKA-a-R-V','PKC-a-b-II_pT638_T641-R-V','PKC-b-II_pS660-R-V','PKC-delta_pS664-R-V','PKCa-R-V','PKM2-R-C','PLC-gamma2_pY759-R-C','PLK1-R-C','PMS2-R-V','PR-R-V','PRAS40_pT246-R-V','PREX1-R-V','PTEN-R-V','PTPN12-R-V','Puma-R-C','PYGB-R-V','Pyk2_pY402-R-C','Rab11-R-C','Rab25-R-V','Rad23A-R-C','Rad50-R-V','Rad51-R-C','Raptor-R-V','RBM15-R-V','Rb_pS807_S811-R-V','Rictor-R-C','Rictor_pT1135-R-V','RIP-R-C','RIP3-R-C','RPA32-R-V','RPA32_pS4_S8-R-C','RRM1-R-C','RRM2-R-C','RSK-R-C','S100A4-R-V','S6_pS235_S236-R-V','S6_pS240_S244-R-V','SDHA-R-V','SETD2-R-C','SFRP1-R-C','Shc_pY317-R-V','SHP-2_pY542-R-C','SHP2-R-V','SLC1A5-R-C','Slfn11-G-C','Smad1-R-V','Smad3-R-V','SOD2-R-V','Sox2-R-V','Src_pY416-R-V','Src_pY527-R-V','Stat3-R-C','Stat3_pY705-R-V','Stat5a-R-V','Stathmin-1-R-V','STING-R-V','TAZ-R-V','TFAM-R-V','TFRC-R-V','TIGAR-R-V','TRIM25-R-C','TRIP13-R-V','TSC1-R-C','TTF1-R-V','Tuberin-R-V','Tuberin_pT1462-R-V','TUFM-R-V','Tyro3-R-V','U-Histone-H2B-R-C','UBAC1-R-V','ULK1_pS757-R-C','UVRAG-R-C','VASP-R-V','VAV1-R-C','VEGFR-2-R-V','Wee1-R-C','Wee1_pS642-R-C','WIPI1-R-C','WIPI2-R-C','XBP-1-G-C','XIAP-R-C','XPF-R-C','XRCC1-R-C','YAP-R-C','YAP_pS127-R-V','YB1_pS102-R-V','ZAP-70-R-C','AceCS1-R-V','ACLY_pS455-R-V','ATM_pS1981-R-V','ATP5H-R-V','ATR-R-V','Aurora-ABC_pT288_pT232_pT198-R-C','Bcl2-R-C','Calnexin-R-V','CD134-R-V','CDK9-R-V','CENP-A-R-V','CHD1L-R-V','CRABP1-R-C','CRABP2-R-V','Cyclin-E1-R-V','DDR1_pY513-R-C','Dvl3-R-V','DYRK1B-R-C','E2F1-R-V','EphA2-R-V','EphA2_pS897-R-C','EphA2_pY588-R-C','ERCC1-R-C','EVI1-R-V','FABP5-R-C','FANCD2-R-V','FRS2-alpha_pY196-R-V','GCLM-R-V','Glucocorticoid-Receptor-R-V','GRB2-R-V','H2AX_pS139-R-C','HER2_pY1248-R-V','Hif-1-alpha-R-C','Histone-H3_pS10-R-V','HLA-DQA1-R-V','HLA-DR-DP-DQ-DX-R-C','HNRNPK-R-V','INPP4b-R-C','IRF-1-R-V','IRF-3-R-V','JNK2-R-V','JNK_pT183_Y185-R-C','KEAP1-R-V','LCN2-R-V','Lyn-R-V','MACC1-R-V','MAPK_pT202_Y204-R-C','MEK1_pS217_S221-R-V','Melanoma-gp100-R-C','Midkine-R-V','MIG6-R-V','NDUFB4-R-V','p27-Kip1-R-C','p38-MAPK-_pT180_Y182-R-V','PAK_pS474_S602_S560-R-V','PAK_pT423_T402-R-C','PAX6-R-V','PAX8-R-C','PERK-R-V','PGM1-R-V','PI3K-p85-R-C','PLC-gamma1-R-V','PLC-gamma1_pS1248-R-V','PLC-gamma2_pY759-R-V','PRC1_pT481-R-C','Rad17_pS645-R-V','RSK1-R-V','SFRP1-R-V','SGK1-R-V','SGK3-R-V','Shc_pY317-R-C','Smad4-R-V','Sox17-R-V','SOX7-R-V','Stat1_pY701-R-V','Stat3_pY705-R-C','Synaptophysin-R-C','TEAD-R-C','VEGFR-2_pY1175-R-C','VHL-R-C','YES1-R-V','ZAP-70-R-V','ZEB1-R-V','ALKBH5-R-V','cGAS-R-V','Cox2-R-V','FTO-R-V','Gli3-R-V','IMP3-R-C','METTL3-R-V','PIP4K2A-R-V','PIP4K2B-R-V','SF2-R-V','SIRP-alpha-R-V','TRIM24-R-C','WTAP-R-V','YTHDF2-R-V','YTHDF3-R-C','Folliculin-R-V','PRMT1-R-V','PRMT5-R-V','Rab11FIP1-R-V','Src-R-V','TSC1-R-V','FGFR1-R-V','FGFR2-R-V','RORA-R-N','TROP2-R-V','Cadherin-6-R-C','Caspase-7-cleaved-R-C','Ezh2-R-V','Folate-Binding-Protein-R-V','VEGFR2_pY1175-R-C')
rppa_to_gene_df <- data.frame(RPPA, Gene)
rppa_to_gene_dict <- setNames(Gene, RPPA)

#Create list of weights for each variant by gene and factor harmonized across all omics
results <- list()
#threshold <- 0.0005
threshold <- -1
get_related_rppa <- function(gene, rppa_to_gene_dict, rppa_df) {
  related_rppa <- names(rppa_to_gene_dict[rppa_to_gene_dict == gene])
  related_rppa <- related_rppa[related_rppa %in% rownames(rppa_df)]
  return(related_rppa)
}

for (factor in colnames(wes_df)) {
  for (variant in rownames(wes_df)) {
    gene <- unlist(strsplit(variant, "_"))[1]
    
    variant_weight <- wes_df[variant, factor]
    gene_weight <- 0
    rppa_weights <- numeric(0)
    
    if (gene %in% rownames(rna_df)) {
      gene_weight <- rna_df[gene, factor]
    }
    
    if (gene %in% rppa_to_gene_dict) {
      related_rppa <- get_related_rppa(gene, rppa_to_gene_dict, rppa_df)
      if (length(related_rppa) > 0) {
        rppa_weights <- sapply(related_rppa, function(rppa) rppa_df[rppa, factor])
      }
    }
    
    #Check GOF conditions
    if ((abs(variant_weight) > threshold & abs(gene_weight) > threshold & 
         sign(variant_weight) == sign(gene_weight)) ||
        (any(abs(variant_weight) > threshold & abs(rppa_weights) > threshold & sign(variant_weight) == sign(rppa_weights)))) {
      results[[length(results) + 1]] <- list(
        Variant = variant,
        Label = "GOF",
        Factor = factor,
        Gene = gene,
        Weights = list(Variant = variant_weight, Gene = gene_weight, RPPA = rppa_weights)
      )
    }
    
    #Check LOF conditions
    if ((abs(variant_weight) > threshold & abs(gene_weight) > threshold & 
         sign(variant_weight) != sign(gene_weight)) ||
        (any(abs(variant_weight) > threshold & abs(rppa_weights) > threshold & sign(variant_weight) != sign(rppa_weights)))) {
      results[[length(results) + 1]] <- list(
        Variant = variant,
        Label = "LOF",
        Factor = factor,
        Gene = gene,
        Weights = list(Variant = variant_weight, Gene = gene_weight, RPPA = rppa_weights)
      )
    }
  }
}


output_df <- data.frame(Variant = character(), Label = character(), Factor = character(), 
                        Variant_Weight = numeric(), Gene_Weight = numeric(), RPPA_Weights = character(), 
                        stringsAsFactors = FALSE)

for (res in results) {
  rppa_weights_str <- if (length(res$Weights$RPPA) > 0) paste(res$Weights$RPPA, collapse = ";") else ""
  
  output_df <- rbind(output_df, data.frame(
    Variant = res$Variant,
    Label = res$Label,
    Factor = res$Factor,
    Variant_Weight = res$Weights$Variant,
    Gene_Weight = res$Weights$Gene,
    RPPA_Weights = rppa_weights_str,
    stringsAsFactors = FALSE
  ))
}

write.csv(output_df, "gof_lof_variants_weights_scaled_all.csv", row.names = FALSE)