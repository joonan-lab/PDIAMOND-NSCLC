#### Code title: run_CNVsignature
#### Description: CNV signature analysis using Sigminer
#### Last log: 23.07.25

## Settings
rm(list = ls())

## Needed libraries
library(tidyverse)
library(NMF)
library(dplyr)
library(sigminer)

# Load CNV data
# Input data should contain six columns: Chromosome, Start.bp, End.bp, modal_cn, minor_cn, and sample
m.cnv = readRDS('files/cnv_facets_merged.Rdata')

# Load sample info and set sex info
s = openxlsx::read.xlsx('files/TableS1_Data_overview.xlsx', sheet = 2)
sex_input = s %>% dplyr::select(c(Sample.ID, Sex))
colnames(sex_input) = c('sample', 'sex')
options(sigminer.sex = sex_input)

# Tally variation components
CNV.obj.combined <- read_copynumber(m.cnv,
                                    genome_build = "hg38",
                                    complement = FALSE, verbose = TRUE,
                                    add_loh = T
)
ncores <- 8

CNV.combined.tally.S <- sig_tally(CNV.obj.combined, method = "S", cores = ncores)

# Estimate the number of signature
EST.combined.S <- sig_estimate(CNV.combined.tally.S$nmf_matrix,
                               range = 2:12, nrun = 50, cores = ncores, use_random = TRUE,
                               save_plots = FALSE,
                               verbose = TRUE
)

p = show_sig_number_survey(EST.combined.S, right_y = NULL)
p = add_h_arrow(p, x = 7.3, y = 0.989, seg_len = 0.5, space = 0.2)
p

# Extract signature profile
cn_sig = sig_extract(CNV.combined.tally.S$all_matrices$CN_48, n_sig = 7)
cosmic_cn = get_sig_similarity(cn_sig, sig_db = "CNS_TCGA")

show_sig_profile(cn_sig,
                 mode = "copynumber",
                 style = "cosmic", method = "S", normalize = "feature")

# Group samples and allocate to each signatures
CNV_group_combined <- get_groups(cn_sig, method = "consensus", match_consensus = TRUE) # The 'enrich_sig' column is set to dominant signature in one group
colnames(CNV_group_combined) <- c("sample", "cnv_group", "cnv_weight", "cnv_enrich_sig")








