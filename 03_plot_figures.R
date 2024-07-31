# Script to look at associations betwen pathways of interest
library(dplyr)
library(ggplot2)
library(corrplot)

selected_pathways <-
  read.csv("/Volumes/auslander/linux/McKenna/pdx/data/mel_pathway/pathways_plot.csv")
selected_pathways <- selected_pathways$pathway

reorder_list <- c("GO_BIOLOGICAL_ADHESION_wes",
                  "GO_BIOLOGICAL_ADHESION_rnaseq",
                  "GO_BIOLOGICAL_ADHESION_rppa",
                  "GO_METAL_ION_TRANSMEMBRANE_TRANSPORTER_ACTIVITY_wes",
                  "GO_METAL_ION_TRANSMEMBRANE_TRANSPORTER_ACTIVITY_rnaseq",
                  "GO_METAL_ION_TRANSMEMBRANE_TRANSPORTER_ACTIVITY_rppa",
                  "GO_NEUROLOGICAL_SYSTEM_PROCESS_wes",
                  "GO_NEUROLOGICAL_SYSTEM_PROCESS_rnaseq",
                  "GO_NEUROLOGICAL_SYSTEM_PROCESS_rppa",
                  "GO_HAIR_CELL_DIFFERENTIATION_wes",
                  "GO_HAIR_CELL_DIFFERENTIATION_rnaseq",
                  "GO_HAIR_CELL_DIFFERENTIATION_rppa",
                  "REACTOME_SIGNALING_BY_VEGF_wes",
                  "REACTOME_SIGNALING_BY_VEGF_rnaseq",
                  "REACTOME_SIGNALING_BY_VEGF_rppa",
                  "REACTOME_DEUBIQUITINATION_wes",
                  "REACTOME_DEUBIQUITINATION_rnaseq",
                  "REACTOME_DEUBIQUITINATION_rppa")

##############################################
############## responders ####################
##############################################
# load in the omics data for responders and nonresponders
responders_corr <-
  read.csv("/wistar/auslander/McKenna/pdx/outputs/mel_pathway/" +
             "for_paper/responders_corr.csv")
rownames(responders_corr) <- responders_corr$X
responders_corr <- responders_corr %>% select(-X)

# filter to only include pathways of interest
pattern <- paste(selected_pathways, collapse = "|")
cols_to_filter <- colnames(responders_corr)
filtered_cols <- cols_to_filter[grepl(pattern, cols_to_filter)]
responders <- responders_corr %>%
  select(all_of(filtered_cols)) # filter columns
# filter rows
responders <- responders[row.names(responders) %in% filtered_cols, ]
responders <- as.matrix(responders)

responders <- responders[rev(reorder_list), rev(reorder_list)]

# make plots
pdf("/wistar/auslander/McKenna/pdx/outputs/mel_pathway/" +
      "wes_rnaseq_rppa/r_corrplot_responders.pdf", width = 10, height = 7)
plot.new()
responders_plot <- corrplot(responders, order = "original",
                            type = "lower", diag = FALSE, tl.cex = 0.5)
dev.off()

##############################################
############## nonresponders #################
##############################################
nonresponders_corr <-
  read.csv("/wistar/auslander/McKenna/pdx/outputs" +
             "/mel_pathway/for_paper/nonresponders_corr.csv")
rownames(nonresponders_corr) <- nonresponders_corr$X
nonresponders_corr <- nonresponders_corr %>% select(-X)

# filter to only include pathways of interest
cols_to_filter <- colnames(nonresponders_corr)
filtered_cols <- cols_to_filter[grepl(pattern, cols_to_filter)]
nonresponders <- nonresponders_corr %>%
  select(all_of(filtered_cols)) # filter columns

# filter rows
nonresponders <- nonresponders[row.names(nonresponders) %in% filtered_cols, ]
nonresponders <- as.matrix(nonresponders)
nonresponders <- nonresponders[rev(reorder_list), rev(reorder_list)]

# make plots
pdf("/wistar/auslander/McKenna/pdx/outputs/mel_pathway/" +
      "wes_rnaseq_rppa/r_corrplot_nonresponders.pdf", width = 10, height = 7)
plot.new()
nonresponders_plot <- corrplot(nonresponders, order = "original",
                               type = "lower", diag = FALSE, tl.cex = 0.5)
dev.off()

##############################################
######### MR/SD grouped together #############
##############################################
mrsd_corr <- read.csv("/wistar/auslander/McKenna/pdx/outputs/" +
                        "mel_pathway/for_paper/mrsd_corr.csv")
rownames(mrsd_corr) <- mrsd_corr$X
mrsd_corr <- mrsd_corr %>% select(-X)

cols_to_filter <- colnames(mrsd_corr)
filtered_cols <- cols_to_filter[grepl(pattern, cols_to_filter)]
mrsd <- mrsd_corr %>%
  select(all_of(filtered_cols)) # filter columns

# filter rows
mrsd <- mrsd[row.names(mrsd) %in% filtered_cols, ]
mrsd <- as.matrix(mrsd)
mrsd <- mrsd[rev(reorder_list), rev(reorder_list)]

# make plots
pdf("/wistar/auslander/McKenna/pdx/outputs/mel_pathway/" +
      "wes_rnaseq_rppa/r_corrplot_mrsd.pdf", width = 10, height = 7)
plot.new()
mrsd_plot <- corrplot(mrsd, order = "original",
                      type = "lower", diag = FALSE, tl.cex = 0.5)
dev.off()


##############################################
############## TT: responders ################
##############################################
# load in the omics data for responders and nonresponders
responders_tt_corr <-
  read.csv("/wistar/auslander/McKenna/pdx/outputs/" +
             "mel_pathway/for_paper/responders_tt_corr.csv")
rownames(responders_tt_corr) <- responders_tt_corr$X
responders_tt_corr <- responders_tt_corr %>% select(-X)

# filter to only include pathways of interest
pattern <- paste(selected_pathways, collapse = "|")
cols_to_filter <- colnames(responders_tt_corr)
filtered_cols <- cols_to_filter[grepl(pattern, cols_to_filter)]
responders_tt <- responders_tt_corr %>%
  select(all_of(filtered_cols)) # filter columns
# filter rows
responders_tt <- responders_tt[row.names(responders_tt) %in% filtered_cols, ]
responders_tt <- as.matrix(responders_tt)

responders_tt <- responders_tt[rev(reorder_list), rev(reorder_list)]

# make plots
pdf("/wistar/auslander/McKenna/pdx/outputs/mel_pathway/" +
      "wes_rnaseq_rppa/r_corrplot_responders_tt.pdf", width = 10, height = 7)
plot.new()
responders_tt_plot <- corrplot(responders_tt, order = "original",
                               type = "lower", diag = FALSE, tl.cex = 0.5)
dev.off()

##############################################
############# TT: nonresponders ##############
##############################################
nonresponders_tt_corr <- read.csv("/wistar/auslander/McKenna/pdx/" +
                                    "outputs/mel_pathway/for_paper/" +
                                    "nonresponders_tt_corr.csv")
rownames(nonresponders_tt_corr) <- nonresponders_tt_corr$X
nonresponders_tt_corr <- nonresponders_tt_corr %>% select(-X)

# filter to only include pathways of interest
cols_to_filter <- colnames(nonresponders_tt_corr)
filtered_cols <- cols_to_filter[grepl(pattern, cols_to_filter)]
nonresponders_tt <- nonresponders_tt_corr %>%
  select(all_of(filtered_cols)) # filter columns

# filter rows
nonresponders_tt <-
  nonresponders_tt[row.names(nonresponders_tt) %in% filtered_cols, ]
nonresponders_tt <- as.matrix(nonresponders_tt)
nonresponders_tt <- nonresponders_tt[rev(reorder_list), rev(reorder_list)]

# make plots
pdf("/wistar/auslander/McKenna/pdx/outputs/mel_pathway/for_paper/" +
      "r_corrplot_nonresponders_tt.pdf", width = 10, height = 7)
plot.new()
nonresponders_tt_plot <- corrplot(nonresponders_tt, order = "original",
                                  type = "lower", diag = FALSE, tl.cex = 0.5)
dev.off()

##############################################
######### TT: MR/SD grouped together #########
##############################################
mrsd_tt_corr <- read.csv("/wistar/auslander/McKenna/pdx/" +
                           "outputs/mel_pathway/for_paper/" +
                           "mrsd_tt_corr.csv")
rownames(mrsd_tt_corr) <- mrsd_tt_corr$X
mrsd_tt_corr <- mrsd_tt_corr %>% select(-X)

cols_to_filter <- colnames(mrsd_tt_corr)
filtered_cols <- cols_to_filter[grepl(pattern, cols_to_filter)]
mrsd_tt <- mrsd_tt_corr %>%
  select(all_of(filtered_cols)) # filter columns

# filter rows
mrsd_tt <- mrsd_tt[row.names(mrsd_tt) %in% filtered_cols, ]
mrsd_tt <- as.matrix(mrsd_tt)
mrsd_tt <- mrsd_tt[rev(reorder_list), rev(reorder_list)]

# make plots
pdf("/wistar/auslander/McKenna/pdx/outputs/mel_pathway/for_paper/" +
      "r_corrplot_mrsd_tt.pdf", width = 10, height = 7)
plot.new()
mrsd_tt_plot <- corrplot(mrsd_tt, order = "original",
                         type = "lower", diag = FALSE, tl.cex = 0.5)
dev.off()


##############################################
#####      RPPA averaged by isoform    #######
##############################################
nonresponders_corr <-
  read.csv("/Volumes/auslander/linux/McKenna/pdx/outputs/mel_pathway/for_paper/isononresponders_corr.csv")
rownames(nonresponders_corr) <- nonresponders_corr$X
nonresponders_corr <- nonresponders_corr %>% select(-X)

# filter to only include pathways of interest
cols_to_filter <- colnames(nonresponders_corr)
filtered_cols <- cols_to_filter[grepl(pattern, cols_to_filter)]
nonresponders <- nonresponders_corr %>%
  select(all_of(filtered_cols)) # filter columns

# filter rows
nonresponders <- nonresponders[row.names(nonresponders) %in% filtered_cols, ]
nonresponders <- as.matrix(nonresponders)
nonresponders <- nonresponders[rev(reorder_list), rev(reorder_list)]

# make plots
pdf("/Volumes/auslander/linux/McKenna/pdx/outputs/mel_pathway/wes_rnaseq_rppa/r_corrplot_isononresponders.pdf", width = 10, height = 7)
nonresponders_plot <- corrplot(nonresponders, order = "original",
                               type = "lower", diag = FALSE, tl.cex = 0.5)
dev.off()

###############################################################################
# responders
responders_corr <-
  read.csv("/Volumes/auslander/linux/McKenna/pdx/outputs/mel_pathway/for_paper/isoresponders_corr.csv")
rownames(responders_corr) <- responders_corr$X
responders_corr <- responders_corr %>% select(-X)

# filter to only include pathways of interest
cols_to_filter <- colnames(responders_corr)
filtered_cols <- cols_to_filter[grepl(pattern, cols_to_filter)]
responders <- responders_corr %>%
  select(all_of(filtered_cols)) # filter columns

# filter rows
responders <- responders[row.names(responders) %in% filtered_cols, ]
responders <- as.matrix(responders)
responders <- responders[rev(reorder_list), rev(reorder_list)]

# make plots
pdf("/Volumes/auslander/linux/McKenna/pdx/outputs/mel_pathway/wes_rnaseq_rppa/r_corrplot_isoresponders.pdf", width = 10, height = 7)
responders_plot <- corrplot(responders, order = "original",
                            type = "lower", diag = FALSE, tl.cex = 0.5)
dev.off()

###############################################################################
# mixed response/stable disease
mrsd_corr <-
  read.csv("/Volumes/auslander/linux/McKenna/pdx/outputs/mel_pathway/for_paper/isomrsd_corr.csv")
rownames(mrsd_corr) <- mrsd_corr$X
mrsd_corr <- mrsd_corr %>% select(-X)

# filter to only include pathways of interest
cols_to_filter <- colnames(mrsd_corr)
filtered_cols <- cols_to_filter[grepl(pattern, cols_to_filter)]
mrsd <- mrsd_corr %>%
  select(all_of(filtered_cols)) # filter columns

# filter rows
mrsd <- mrsd[row.names(mrsd) %in% filtered_cols, ]
mrsd <- as.matrix(mrsd)
mrsd <- mrsd[rev(reorder_list), rev(reorder_list)]

# make plots
pdf("/Volumes/auslander/linux/McKenna/pdx/outputs/mel_pathway/wes_rnaseq_rppa/r_corrplot_isomrsd.pdf", width = 10, height = 7)
mrsd_plot <- corrplot(mrsd, order = "original",
                      type = "lower", diag = FALSE, tl.cex = 0.5)
dev.off()