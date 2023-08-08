#### Code title: run_viper
#### Description: Inferred protein activity based on the transcriptome data
#### Last log: 23.08.06

## Needed libraries
library(viper)
library(dplyr)
library(plyr)
library(aracne.networks)
library(Biobase)
library(bigreadr)
library(stringr)
library(data.table)
library('org.Hs.eg.db')

# Load data
data(regulonluad)
data(regulonlusc)
data(package="aracne.networks")$result[,"Item"]

ad_data <- read.csv("AD_expression_data.csv", row.names = 1, header=T)
sc_data <- read.csv("SC_expression_data.csv", row.names = 1, header=T)
satpathy_data <- read.csv("satpathy_expression_data.csv", row.names = 1, header=T)
gillette_data <- read.csv("gillette_expression_data.csv", row.names = 1, header=T)
# "AD, SC_expression_data.csv"; Row: gene symbol x Column: Sample, RNA dataframe, Input data represent RNA expression values of patients in our cohort.
# "satpathy, gillette_expression_data.csv"; Row: gene symbol x Column: Sample, RNA dataframe, Input data correspond to RNA expression values from "Satpathy et al" and "Gillette et al," respectively.


# gene map function
entrez_id <- keys(org.Hs.eg.db, keytype = "ENTREZID")
org <- select(org.Hs.eg.db, keys = entrez_id, keytype = "ENTREZID", columns = "SYMBOL")
colnames(org) <- c("gene_id", "symbol")

# id → symbol
converter <- function(st){
  return(org[org$gene_id==st,]$symbol)
}

# symbol → id
s_2_i <- function(st){
  return(org_f[org_f$symbol==st,]$gene_id)
}

his <- "AD"
data1_gene <- rownames(ad_data)
org_gene <- c(org$symbol)

data1_f <- ad_data[intersect(data1_gene, org_gene),]
org_f <- org[org$symbol %in% intersect(data1_gene, org_gene),]
data1_gene <- rownames(data1_f)
length(data1_gene)

rownames(data1_f) <- sapply(data1_gene, s_2_i)

vpres = viper(data1_f, regulonluad, verbose = FALSE) ## single sample viper
geneid_ss = ldply(sapply(row.names(vpres), converter), data.frame) ## gene mapping
colnames(geneid_ss) = c('id', 'gene')
vpres_f <- vpres[geneid_ss$id,]
rownames(vpres_f) <- geneid_ss$gene

his <- "SC"
sc_data[1:5,1:5]
data1_gene <- rownames(sc_data)
org_gene <- c(org$symbol)

data1_f <- sc_data[intersect(data1_gene, org_gene),]
org_f <- org[org$symbol %in% intersect(data1_gene, org_gene),]
data1_gene <- rownames(data1_f)

rownames(data1_f) <- sapply(data1_gene, s_2_i)

vpres = viper(data1_f, regulonlusc, verbose = FALSE) ## single sample viper
geneid_ss = ldply(sapply(row.names(vpres), converter), data.frame) ## gene mapping
colnames(geneid_ss) = c('id', 'gene')
vpres_f <- vpres[geneid_ss$id,]
rownames(vpres_f) <- geneid_ss$gene


## CPTAC LUAD and LSCC datasets.
his <- "satpathy"
rownames(satpathy_data) <- satpathy_data[,1]
satpathy_data <- satpathy_data[,-1]

data1_gene <- rownames(satpathy_data)
org_gene <- c(org$symbol)

data1_f <- satpathy_data[intersect(data1_gene, org_gene),]
org_f <- org[org$symbol %in% intersect(data1_gene, org_gene),]
data1_gene <- rownames(data1_f)

rownames(data1_f) <- sapply(data1_gene, s_2_i)

vpres = viper(data1_f, regulonlusc, verbose = FALSE) ## single sample viper
geneid_ss0 = ldply(sapply(row.names(vpres), converter), data.frame) ## gene mapping
colnames(geneid_ss0) = c('id', 'gene')
vpres_f <- vpres[geneid_ss0$id,]
rownames(vpres_f) <- geneid_ss0$gene


his <- "gillette"
rownames(gillette_data) <- gillette_data[,1]
gillette_data <- gillette_data[,-1]

data1_gene <- rownames(gillette_data)
org_gene <- c(org$symbol)

data1_f <- gillette_data[intersect(data1_gene, org_gene),]
org_f <- org[org$symbol %in% intersect(data1_gene, org_gene),]
data1_gene <- rownames(data1_f)

rownames(data1_f) <- sapply(data1_gene, s_2_i)

vpres = viper(data1_f, regulonluad, verbose = FALSE) ## single sample viper
geneid_ss0 = ldply(sapply(row.names(vpres), converter), data.frame) ## gene mapping
colnames(geneid_ss0) = c('id', 'gene')
vpres_f <- vpres[geneid_ss0$id,]
rownames(vpres_f) <- geneid_ss0$gene