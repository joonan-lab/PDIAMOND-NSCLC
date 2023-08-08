#### Code title: run_regression
#### Description: Linear regression between immune cells or the status of immune clusters and gene expression
#### Last log: 23.08.06

## Needed libraries
library(dplyr)
library(progress)
library(data.table)

# Load data
rna_data <- read.csv("rna_expression_data.csv", header=T)
protein_exp_data <- read.csv("protein_expression_data.csv", header=T)
protein_activity_data <- read.csv("protein_activity_data.csv", header=T)
# Three dataset formats; Row: Samples x Column: 11 major immune or stromal cells and RNA expression values


# run regression from transcriptome data
features = colnames(rna_data[,2:11])
genes = colnames(rna_data[,21:ncol(rna_data)])

heatmap_df = c()
p_df = c()
fdr_df = c()
filtered_genes = c()

pb <- progress_bar$new(
  format = " Progress: [:bar] :percent, Estimated completion time: :eta",
  total = 10000, # totalnumber of ticks to complete (default 100)
  clear = FALSE, # whether to clear the progress bar on completion (default TRUE)
  width= 80 # width of the progress bar
)
for (i in features){
  pb$tick()
  coeffs = c()
  pvals = c()
  for (gene in genes){
    folm = paste0(i, "~",gene,"+DX+rna_batch")
    result = lm(folm, data = rna_data)
    coeffs = c(coeffs, summary(result)$coefficient[2])
    pvals = c(pvals, summary(result)$coefficient[14])
  }
  heatmap_df <- cbind(heatmap_df, coeffs)
  p_df <- cbind(p_df, pvals)
  fdrs <- p.adjust(pvals, "fdr")
  fdr_df <- cbind(fdr_df, fdrs)
  res <- data.frame(genes, coeffs, pvals)
  res$fdr <- p.adjust(res$pvals, "fdr")
  res$scaled_fdr <- -log10(res$fdr)
  res$scaled_pvals <- -log10(res$pvals)
}
colnames(heatmap_df) <- features
rownames(heatmap_df) <- genes

colnames(p_df) <- features
rownames(p_df) <- genes

colnames(fdr_df) <- features
rownames(fdr_df) <- genes

heatmap_df <- as.data.frame(heatmap_df)
p_df <- as.data.frame(p_df)
fdr_df <- as.data.frame(fdr_df)

coeffs = c()
pvals = c()
for (gene in genes){
  folm = paste0("xcell_subtype~",gene,"+DX+rna_batch")
  result = glm(folm, data=rna_data, family = binomial)
  coeffs = c(coeffs, summary(result)$coefficient[2])
  pvals = c(pvals, summary(result)$coefficient[14])
}
fdrs = p.adjust(pvals, "fdr")
res <- data.frame(genes, coeffs, pvals)

heatmap_df$xcell_subtype <- coeffs
p_df$xcell_subtype <- pvals
fdr_df$xcell_subtype <- fdrs

heat_genes = c()
masked_df <- data.frame(heatmap_df)
fdr_genes = c()

# p_value > 0.05 → NA
for (i in colnames(p_df)){
  heat_genes <- rownames(subset(p_df, p_df[,i] > 0.05))
  masked_df[heat_genes, i] = NA
}

filtered_df = subset(masked_df, rowSums(abs(masked_df), na.rm = T)!=0)
filtered_df = filtered_df[,colSums(is.na(filtered_df))<nrow(filtered_df)]

pheat_df <- data.frame(filtered_df)
pheat_df[is.na(pheat_df)] <- 0


# Protein exp
features = colnames(protein_exp_data[,2:11])
genes = colnames(protein_exp_data[,21:ncol(protein_exp_data)])

heatmap_df = c()
p_df = c()
fdr_df = c()
filtered_genes = c()

pb <- progress_bar$new(
  format = " Progress: [:bar] :percent, Estimated completion time: :eta",
  total = 10000, # totalnumber of ticks to complete (default 100)
  clear = FALSE, # whether to clear the progress bar on completion (default TRUE)
  width= 80 # width of the progress bar
)
for (i in features){
  pb$tick()
  coeffs = c()
  pvals = c()
  for (gene in genes){
    folm = paste0(i, "~",gene,"+DX+rna_batch")
    result = lm(folm, data = protein_exp_data)
    coeffs = c(coeffs, summary(result)$coefficient[2])
    pvals = c(pvals, summary(result)$coefficient[14])
  }
  heatmap_df <- cbind(heatmap_df, coeffs)
  p_df <- cbind(p_df, pvals)
  fdrs <- p.adjust(pvals, "fdr")
  fdr_df <- cbind(fdr_df, fdrs)
  res <- data.frame(genes, coeffs, pvals)
  res$fdr <- p.adjust(res$pvals, "fdr")
  res$scaled_fdr <- -log10(res$fdr)
  res$scaled_pvals <- -log10(res$pvals)
}
colnames(heatmap_df) <- features
rownames(heatmap_df) <- genes

colnames(p_df) <- features
rownames(p_df) <- genes

colnames(fdr_df) <- features
rownames(fdr_df) <- genes

heatmap_df <- as.data.frame(heatmap_df)
p_df <- as.data.frame(p_df)
fdr_df <- as.data.frame(fdr_df)

coeffs = c()
pvals = c()
for (gene in genes){
  folm = paste0("xcell_subtype~",gene,"+DX+rna_batch")
  result = glm(folm, data=protein_exp_data, family = binomial)
  coeffs = c(coeffs, summary(result)$coefficient[2])
  pvals = c(pvals, summary(result)$coefficient[14])
}
fdrs = p.adjust(pvals, "fdr")
res <- data.frame(genes, coeffs, pvals)
write.csv(res, 'xcell_protein_exp_regression_dx+batch.csv')

heatmap_df$xcell_subtype <- coeffs
p_df$xcell_subtype <- pvals
fdr_df$xcell_subtype <- fdrs

heat_genes = c()
masked_df <- data.frame(heatmap_df)
fdr_genes = c()

# p_value > 0.05 → NA
for (i in colnames(p_df)){
  heat_genes <- rownames(subset(p_df, p_df[,i] > 0.05))
  masked_df[heat_genes, i] = NA
}

filtered_df = subset(masked_df, rowSums(abs(masked_df), na.rm = T)!=0)
filtered_df = filtered_df[,colSums(is.na(filtered_df))<nrow(filtered_df)]

pheat_df <- data.frame(filtered_df)
pheat_df[is.na(pheat_df)] <- 0


# Protein activity
features = colnames(protein_act_data[,2:11])
genes = colnames(protein_act_data[,21:ncol(protein_act_data)])

heatmap_df = c()
p_df = c()
fdr_df = c()
filtered_genes = c()

pb <- progress_bar$new(
  format = " Progress: [:bar] :percent, Estimated completion time: :eta",
  total = 10000, # totalnumber of ticks to complete (default 100)
  clear = FALSE, # whether to clear the progress bar on completion (default TRUE)
  width= 80 # width of the progress bar
)
for (i in features){
  pb$tick()
  coeffs = c()
  pvals = c()
  for (gene in genes){
    folm = paste0(i, "~",gene,"+DX+rna_batch")
    result = lm(folm, data = protein_act_data)
    coeffs = c(coeffs, summary(result)$coefficient[2])
    pvals = c(pvals, summary(result)$coefficient[14])
  }
  heatmap_df <- cbind(heatmap_df, coeffs)
  p_df <- cbind(p_df, pvals)
  fdrs <- p.adjust(pvals, "fdr")
  fdr_df <- cbind(fdr_df, fdrs)
  res <- data.frame(genes, coeffs, pvals)
  res$fdr <- p.adjust(res$pvals, "fdr")
  res$scaled_fdr <- -log10(res$fdr)
  res$scaled_pvals <- -log10(res$pvals)
}
colnames(heatmap_df) <- features
rownames(heatmap_df) <- genes

colnames(p_df) <- features
rownames(p_df) <- genes

colnames(fdr_df) <- features
rownames(fdr_df) <- genes

heatmap_df <- as.data.frame(heatmap_df)
p_df <- as.data.frame(p_df)
fdr_df <- as.data.frame(fdr_df)

coeffs = c()
pvals = c()
for (gene in genes){
  folm = paste0("xcell_subtype~",gene,"+DX+rna_batch")
  result = glm(folm, data=protein_act_data, family = binomial)
  coeffs = c(coeffs, summary(result)$coefficient[2])
  pvals = c(pvals, summary(result)$coefficient[14])
}
fdrs = p.adjust(pvals, "fdr")
res <- data.frame(genes, coeffs, pvals)

heatmap_df$xcell_subtype <- coeffs
p_df$xcell_subtype <- pvals
fdr_df$xcell_subtype <- fdrs

heat_genes = c()
masked_df <- data.frame(heatmap_df)
fdr_genes = c()

# p_value > 0.05 → NA
for (i in colnames(p_df)){
  heat_genes <- rownames(subset(p_df, p_df[,i] > 0.05))
  masked_df[heat_genes, i] = NA
}

filtered_df = subset(masked_df, rowSums(abs(masked_df), na.rm = T)!=0)
filtered_df = filtered_df[,colSums(is.na(filtered_df))<nrow(filtered_df)]

pheat_df <- data.frame(filtered_df)
pheat_df[is.na(pheat_df)] <- 0