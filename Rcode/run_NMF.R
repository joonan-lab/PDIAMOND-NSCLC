#### Code title: run_NMF
#### Description: Make input matrix for NMF using our study multi-omics data
#### Last log: 23.07.14

## Settings
rm(list = ls())

## Needed libraries
library(tidyverse)
library(impute)
library(matrixStats)
library(NMF)


## Load Proteomics data (You can downloaded preprocessed proteomic data from supplementary files)
i_glo = readxl::read_xlsx('TableS7', sheet = 3) %>% column_to_rownames('ID') %>% dplyr::select(contains('-T'))
i_pho = readxl::read_xlsx('TableS7', sheet = 4) %>% column_to_rownames('ID') %>% dplyr::select(contains('-T'))
i_ace = readxl::read_xlsx('TableS7', sheet = 5) %>% column_to_rownames('ID') %>% dplyr::select(contains('-T'))


## Function for make Z-score matrix for NMF
z_NMF <- function(x, missing_cut, neighbors, variance_cut){
  
  ## Imputation
  if (!is.null(missing_cut)){
    print ('Do impute')
    # Remove features with many NAs
    temp1 = x[rowSums(is.na(x)) <= (ncol(x)*missing_cut), ]
    
    # Impute missing expression data
    rn = rownames(temp1)
    temp1_1 = temp1[rowSums(is.na(temp1))!=0,] # for imputation
    temp1_2 = temp1[rowSums(is.na(temp1))==0,]
    
    imp.res = impute.knn(temp1_1 %>% as.matrix(), k = neighbors, rowmax = missing_cut)
    temp1_3 <- imp.res[["data"]] 
    temp1 = rbind.data.frame(temp1_2, temp1_3)
    temp1 = temp1[rn,]
    
  } else {
    print ('No needed to impute')
    temp1 = x %>% as.matrix()
    temp1 = temp1[rowSums(is.na(temp1)) == 0,]
    temp1 = temp1[rowSums(temp1)!=0,]
  }
  
  # Remove low variance
  sd1 <- apply(temp1, 1, sd)
  temp1 = temp1[sd1 > as.numeric(quantile(sd1, variance_cut)),]
  print (paste('Total Features after variance cut', nrow(temp1)))
  
  # Z transformation
  temp1 = base::scale(t(temp1))
  temp1 = as.data.frame(t(temp1))
  colnames(temp1) <- gsub('.', '-', colnames(temp1), fixed = T)
  
  # Create a signed matrix
  temp2 <- temp1
  row.names(temp1) <- paste(row.names(temp1), '1', sep = '_')
  row.names(temp2) <- paste(row.names(temp2), '2', sep = '_')
  
  # 1 : Create one data matrix with all negative numbers zeroed
  temp1[temp1 < 0] <- 0
  # 2 : Create another data matrix with all positive numbers zeroed and
  #     the signs of all negative numbers removed
  temp2[temp2 > 0] <- 0
  temp2 <- temp2*(-1)
  # 3 : Concatenate both matrices resulting in a data matrix twice as large as the original
  res <- rbind(temp1, temp2)
  
  # Remove Null values (all zero)
  res = res[rowSums(res) != 0, ]
  
  return(res)
}

## Run function
z_glo <- z_NMF(i_glo, missing_cut, neighbors, variance_cut)
z_pho <- z_NMF(i_pho, missing_cut, neighbors, variance_cut)
z_ace <- z_NMF(i_ace, missing_cut, neighbors, variance_cut)

## Make NMF input file (GPA:  Global Protein + Phospho + Acetyl)
s_common = intersect(intersect(colnames(d_gl), colnames(d_pp)), colnames(d_ac))
d <- rbind(z_glo[,s_common], z_pho[,s_common], z_ace[,s_common])


## Run NMF
args<-commandArgs(TRUE)
remotes::install_github("renozao/NMF@0.30.4.900")
sessionInfo()

n_run = 200 # If do optimal test, set 50.
mi = 5000 # If do optimal test, set 2000.
n_rank = 5 # Test of optimal rank for data should be first.

if (mi == 'No'){
  # Finding optimal rank
  res <- nmf(d, seq(2,10), method='brunet', nrun=n_run, seed=42, .opt='v2')
  
} else{
  # If you find optimal rank, put your rank above 'Set up the run' part
  res <- nmf(d, n_rank, method='brunet', nrun=n_run, seed=42, .opt='v2', maxIter=as.integer(mi))
}

## Heatmap
s = readxl::read_xlsx('files/TableS1_Data_overview.xlsx', sheet = 2)
annotation_color <- list(
  Sex = c('F' = '#a4e2f4', 'M' = '#740cbd'),
  Smoking_Hx_N_Ex_Cr = c('Cr' = '#800000', 'Ex' = '#FFA500', 'N' = '#E6E6E6', 'NA' = '#BEBEBE'),
  DX = c("AD" = '#b2df8a', 
         "SC" = '#ffa500',
         "NC" = "#feb5da",
         "MA" = "#766aa5",
         "Others" = "#009292"),
  AJCC.8th.TNMstage = c("IA1" = "#fde725" , "IA2" = '#e2e418', "IA3" = '#c5e021', "IB" = '#a8db34', 
                        'IIA' = '#42be71', 'IIB' = '#2fb47c', 
                        'IIIA' = '#26818e', 'IIIB' = '#2f6c8e', 
                        'IVA' = '#440154'),
  TNMstage = c("IA" = "#fde725" , "IB" = '#a8db34', 'IIA' = '#42be71', 'IIB' = '#2fb47c', 'IIIA' = '#26818e', 'IIIB' = '#2f6c8e', 'IVA' = '#440154'),
  TIL.pattern = c('Unknown' = 'grey',
                  'absent' = '#f6ebc8',
                  'non-brisk multifocal' = '#bda881',
                  'brisk band-like' = '#6c4b25',
                  'brisk diffuse' = '#522f09'),
  Pathologic.N = c('N0' = '#E6E6E6', 'N1' = '#9a9945', 'N2' = '#c76674'),
  Recur.status = c('NA' = 'darkgrey', '0' = 'grey', '1' = '#FF7777'),
  Adjuvant.Treatment = c('None' = 'grey', 'CTx' = 'dodgerblue', 'RTx' = 'Firebrick', 'CTx & RTx' = 'purple'),
  Subtype = c('1' = '#fb8072', '2' = '#bebada', '3' = '#ffd61f', '4' = '#8dd3c7', '5' = '#80b1d3'),
  cluster_core = c('Y'='black', 'N'='white'),
  WGD = c('Y'='black', 'N'='white'),
  TP53 = c('Not.altered'='white', 'CNV.loss'='#FA8072', 'Mutation'='#FFDEAD', 'Mut+CNV.loss'='#A0522D'),
  Others_TSG = c('Not.altered'='white', 'CNV.loss'='#FA8072', 'Mutation'='#FFDEAD', 'Mut+CNV.loss'='#A0522D'),
  EGFR = c('None'='white','Others_Indel'='#ff7008','Others_SNV'='#338cc7','L858R'='#b0bf1a', 'exon19_del'='#b01f35'),
  Others_Oncogene = c('ERBB2'='#4e0090', 'KRAS'='#01524a', 'PIK3CA'='#ff6e12','KRAS_PIK3CA'='#a14107','MET'='#5a5a5a',
                      'ALK'='#a12049','ROS1'='#19189e','RET'='#431f16',
                      'None'='white'),
  celltype_based_subtype = c('Cold_Immunogram'= 'blue', 'Hot_Immunogram' ='firebrick', 'NA' = 'grey')
  )

annot.mat = s %>%
  mutate(cluster_core = ifelse(Subtype.membership >= 0.5, 'Y', 'N')) %>%
  dplyr::select(
    'Others_Oncogene'=Other.Oncogene, EGFR,'Others_TSG'= Other.TSG, TP53, WGD,
    'celltype_based_subtype' = Celltype.based.subtype, TIL.pattern,
    Recur.status, Adjuvant.Treatment,
    Pathologic.N, AJCC.8th.TNMstage,
    'Smoking_Hx_N_Ex_Cr' = Smoking.history,Sex, Age, DX, 
    Subtype.membership, cluster_core, Subtype, Sample.ID) %>% column_to_rownames('Sample.ID') %>%
  arrange(Subtype, desc(Subtype.membership))

mat = as.matrix(as.data.frame(res@consensus))
colnames(mat) <- rownames(mat) <- gsub('.', '-', rownames(mat), fixed = T)
mat = mat[rownames(annot.mat),rownames(annot.mat)]

pheatmap::pheatmap(mat, cluster_cols = F, legend=T, cluster_rows = F,
                   border_color = 'white',
                   gaps_col = c(55,100,152,195,229),
                   annotation_col = annot.mat,
                   annotation_colors = annotation_color,
                   show_rownames = F, show_colnames = F,
                   color=colorRampPalette(c("navy", "white", "red"))(20))
