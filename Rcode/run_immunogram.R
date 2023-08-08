#### Code title: run_Immunogram
#### Description: Pathway analysis through immune-related pathways using gene expression data
#### Last log: 23.08.06

## Needed libraries
library(GSVA)
library(dplyr)
library(pathwayPCA)

# Load daata
expr <- read.csv("tpm_expression_data.csv", row.names=1)
img_gene <- read.table("immunogram_gene_list.txt", sep='\t', header=TRUE) # immunogram gene list
# "tpm_expression_data.csv"; Row: gene symbol x Column: Sample, RNA dataframe, Input data should only contain information about tumor samples.
# "immunogram_gene_list.txt"; Input data should contain two columns: gene_sets and gene, Immunogram gene list can be obtained and utilized by accessing the Karasaki et al. "An immunogram for the cancer-immunity cycle: towards personalized immunotherapy of lung cancer."


# Run Immunogram
expr <- as.matrix(expr)
img_gene$Gene..set <- rownames(img_gene)
img_gene <- img_gene[,1:2]

# igs axis 7
t_cell_immunity = c(img_gene[img_gene$Gene..set == "T cell (Immune metagene)", 2])
activated_DC = c(img_gene[img_gene$Gene..set == "activated DC (LM22)", 2])
inhibitor_cells = c(img_gene[img_gene$Gene..set == "inhibitor cells", 2])
trafficking = c(img_gene[img_gene$Gene..set == "Trafficking & infiltration", 2])
recognition = c(img_gene[img_gene$Gene..set == "Recognition of tumor cells", 2])
immune_checkpoints = c(img_gene[img_gene$Gene..set == "Immune checkpoints", 2])
other_inhibitory = c(img_gene[img_gene$Gene..set == "Other inhibitory molecules", 2])

myPathways_ls <- list(
  pathway1 = t_cell_immunity,
  pathway2 = activated_DC,
  pathway3 = inhibitor_cells,
  pathway4 = trafficking,
  pathway5 = recognition,
  pathway6 = immune_checkpoints,
  pathway7 = other_inhibitory
)

myPathway_names <- c(
  "T cell immunity",
  "Priming and activation",
  "Inhibitor cells",
  "Trafficking and infiltration",
  "Recognition of tumor cells",
  "Checkpoint expression",
  "Inhibitory molecules"
)

pathway_list <- CreatePathwayCollection(
  sets_ls = myPathways_ls,
  TERMS = myPathway_names
)

# run GSVA
img_gsva <- gsva(expr, pathway_list$pathways, method="gsva")