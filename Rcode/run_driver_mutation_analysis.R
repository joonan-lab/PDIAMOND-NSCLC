#### Code title: run_linear_regression
#### Description: run linear regression about driver mutation gene
#### Last log: 23.08.06

## Needed libraries
library(data.table)

# Load data
rna_data <- data.frame(fread("RNA_expression_data.csv", header = TRUE))
protein_exp_data <- data.frame(fread("Protein_expression_data.csv"), row.names = 1)
protein_act_data <- data.frame(fread("Protein_activity_data.csv"), row.names = 1)
# Three dataset formats; Row: Samples x Column: status of driver mutation genes and expression values of immunomodulator genes.
# For driver mutations, they are represented in the form of "mutation_x"


find_driver_gene <- function(df){
  
  d_gene <- c()
  for (i in 1:length(colnames(df))){
    if (strsplit(colnames(df)[i], "_")[[1]][1] == "mutation"){
      d_gene <- c(d_gene, colnames(data)[i])
    }
  }
  return (d_gene)
}

row.names(rna_data) <- rna_data[,1]
rna_data <- rna_data[,-1]

d_gene <- find_driver_gene(rna_data) # driver gene
i_gene <- colnames(rna_data[(length(d_gene)+4):ncol(rna_data)]) # immunomodulator gene

# z-score normalization
data_f <- as.data.frame(rna_data)
for (i in i_gene){
  v_l <- data_f[,i]
  v_l_f <- (v_l - mean(v_l)) / sd(v_l)

  data_f[,i] <- v_l_f
}

data_f <- data_f[,i_gene]
data <- cbind(data_f, rna_data[,d_gene])

heatmap_df = c()
p_df = c()
fdr_df = c()
filtered_genes = c()
for (d in d_gene){
  coeffs <- c()
  pvals <- c()
  for (i in i_gene){
    folm <- paste0(i, "~", d)
    result <- lm(folm, data = data)
    
    coeffs <- c(coeffs, summary(result)$coefficient[2])
    pvals <- c(pvals, summary(result)$coefficient[8])
    
    coeff_df <- data.frame(coeffs)
    pval_df <- data.frame(pvals)
  }
  heatmap_df <- cbind(heatmap_df, coeffs)
  p_df <- cbind(p_df, pvals)
  fdrs <- p.adjust(pvals, "fdr")
  fdr_df <- cbind(fdr_df, fdrs)
  res <- data.frame(i_gene, coeffs, pvals)
  res$fdr <- p.adjust(res$pvals, "fdr")
  res$scaled_fdr <- -log10(res$fdr)
  res$scaled_pvals <- -log10(res$pvals)
}
colnames(heatmap_df) <- d_gene
rownames(heatmap_df) <- i_gene

min(heatmap_df)
max(heatmap_df)

colnames(p_df) <- d_gene
rownames(p_df) <- i_gene

colnames(fdr_df) <- d_gene
rownames(fdr_df) <- i_gene

masked_df <- as.data.frame(heatmap_df)
p_df <- as.data.frame(p_df) 
fdr_df <- as.data.frame(fdr_df)
masked_df[1:5,1:5]

# p_value > 0.05 → NA
for (i in colnames(p_df)){
  heat_genes <- rownames(subset(p_df, p_df[,i] > 0.05))
  masked_df[heat_genes, i] = NA
}
masked_df[1:5,1:5]

# fdr > 0.05 → NA
masked_df <- as.data.frame(heatmap_df)
masked_df[1:5, 1:5]
for (i in colnames(fdr_df)){
  heat_genes <- rownames(subset(fdr_df, fdr_df[,i] > 0.05))
  masked_df[heat_genes, i] = NA
}
masked_df[1:5, 1:5]

# fdr > 0.01 → NA
masked_df <- as.data.frame(heatmap_df)
masked_df[1:5, 1:5]
for (i in colnames(fdr_df)){
  heat_genes <- rownames(subset(fdr_df, fdr_df[,i] > 0.01))
  masked_df[heat_genes, i] = NA
}
masked_df[1:5, 1:5]

## protein exp
d_gene <- find_driver_gene(protein_exp_data) # driver gene
i_gene <- colnames(protein_exp_data[(length(d_gene)+4):ncol(protein_exp_data)]) # immunomodulator gene

data_f <- as.data.frame(protein_exp_data)
data_f <- data_f[,i_gene]
data_f[is.na(data_f)] <- 0

for (i in i_gene){
  v_l <- data_f[,i]
  v_l_f <- (v_l - mean(v_l)) / sd(v_l)
  
  data_f[,i] <- v_l_f
}
data_f <- data_f[,i_gene]
data <- cbind(data_f, protein_exp_data[,d_gene])

heatmap_df = c()
p_df = c()
fdr_df = c()
filtered_genes = c()
for (d in d_gene){
  coeffs <- c()
  pvals <- c()
  for (i in i_gene){
    folm <- paste0(i, "~", d)
    result <- lm(folm, data = data)
    
    coeffs <- c(coeffs, summary(result)$coefficient[2])
    pvals <- c(pvals, summary(result)$coefficient[8])
    
    coeff_df <- data.frame(coeffs)
    pval_df <- data.frame(pvals)
  }
  heatmap_df <- cbind(heatmap_df, coeffs)
  p_df <- cbind(p_df, pvals)
  fdrs <- p.adjust(pvals, "fdr")
  fdr_df <- cbind(fdr_df, fdrs)
  res <- data.frame(i_gene, coeffs, pvals)
  res$fdr <- p.adjust(res$pvals, "fdr")
  res$scaled_fdr <- -log10(res$fdr)
  res$scaled_pvals <- -log10(res$pvals)
}
colnames(heatmap_df) <- d_gene
rownames(heatmap_df) <- i_gene

min(heatmap_df)
max(heatmap_df)

colnames(p_df) <- d_gene
rownames(p_df) <- i_gene

colnames(fdr_df) <- d_gene
rownames(fdr_df) <- i_gene

masked_df <- as.data.frame(heatmap_df)
p_df <- as.data.frame(p_df) 
fdr_df <- as.data.frame(fdr_df)

# p_value > 0.05 → NA
for (i in colnames(p_df)){
  heat_genes <- rownames(subset(p_df, p_df[,i] > 0.05))
  masked_df[heat_genes, i] = NA
}

# fdr > 0.05 → NA
masked_df <- as.data.frame(heatmap_df)
masked_df[1:5, 1:5]
for (i in colnames(fdr_df)){
  heat_genes <- rownames(subset(fdr_df, fdr_df[,i] > 0.05))
  masked_df[heat_genes, i] = NA
}

# fdr > 0.01 → NA
masked_df <- as.data.frame(heatmap_df)
masked_df[1:5, 1:5]
for (i in colnames(fdr_df)){
  heat_genes <- rownames(subset(fdr_df, fdr_df[,i] > 0.01))
  masked_df[heat_genes, i] = NA
}

## protein act
d_gene <- find_driver_gene(protein_act_data) # driver gene
i_gene <- colnames(data[(length(d_gene)+4):ncol(protein_act_data)]) # immunomodulator gene

data_f <- as.data.frame(protein_act_data)
data_f <- data_f[,i_gene]

data_f[is.na(data_f)] <- 0

for (i in i_gene){
  v_l <- data_f[,i]
  v_l_f <- (v_l - mean(v_l)) / sd(v_l)
  
  data_f[,i] <- v_l_f
}
data_f <- data_f[,i_gene]
data <- cbind(data_f, protein_act_data[,d_gene])

heatmap_df = c()
p_df = c()
fdr_df = c()
filtered_genes = c()
for (d in d_gene){
  coeffs <- c()
  pvals <- c()
  for (i in i_gene){
    folm <- paste0(i, "~", d)
    result <- lm(folm, data = data)
    
    coeffs <- c(coeffs, summary(result)$coefficient[2])
    pvals <- c(pvals, summary(result)$coefficient[8])
    
    coeff_df <- data.frame(coeffs)
    pval_df <- data.frame(pvals)
  }
  heatmap_df <- cbind(heatmap_df, coeffs)
  p_df <- cbind(p_df, pvals)
  fdrs <- p.adjust(pvals, "fdr")
  fdr_df <- cbind(fdr_df, fdrs)
  res <- data.frame(i_gene, coeffs, pvals)
  res$fdr <- p.adjust(res$pvals, "fdr")
  res$scaled_fdr <- -log10(res$fdr)
  res$scaled_pvals <- -log10(res$pvals)
}
colnames(heatmap_df) <- d_gene
rownames(heatmap_df) <- i_gene

min(heatmap_df)
max(heatmap_df)

colnames(p_df) <- d_gene
rownames(p_df) <- i_gene

colnames(fdr_df) <- d_gene
rownames(fdr_df) <- i_gene

masked_df <- as.data.frame(heatmap_df)
p_df <- as.data.frame(p_df) 
fdr_df <- as.data.frame(fdr_df)

# p_value > 0.05 → NA
for (i in colnames(p_df)){
  heat_genes <- rownames(subset(p_df, p_df[,i] > 0.05))
  masked_df[heat_genes, i] = NA
}

# fdr > 0.05 → NA
masked_df <- as.data.frame(heatmap_df)
masked_df[1:5, 1:5]
for (i in colnames(fdr_df)){
  heat_genes <- rownames(subset(fdr_df, fdr_df[,i] > 0.05))
  masked_df[heat_genes, i] = NA
}

# fdr > 0.01 → NA
masked_df <- as.data.frame(heatmap_df)
masked_df[1:5, 1:5]
for (i in colnames(fdr_df)){
  heat_genes <- rownames(subset(fdr_df, fdr_df[,i] > 0.01))
  masked_df[heat_genes, i] = NA
}