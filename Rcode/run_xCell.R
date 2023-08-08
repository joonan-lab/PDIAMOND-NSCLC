#### Code title: run_xCell
#### Description: Analysis based on cell type enrichment from gene expression data
#### Last log: 23.08.06

## Needed libraries
library(xCell)
library(CancerSubtypes)
library(data.table)

# Load data
tpm_denormalization <- read.csv("tpm_expression_data.csv", row.names=1)
# "tpm_expression_data.csv"; Row: gene symbol x Column: Sample, RNA dataframe
# Input data includes information about both tumor samples and normal samples.


# run xCell
# To obatin raw enrichment score between sample and immune or stroma cell
rawEnrichmentAnalysis <- function(expr, signatures, genes, file.name = NULL, parallel.sz = 4, 
                                  parallel.type = 'SOCK') {
  
  # Reduce the expression dataset to contain only the required genes
  shared.genes <- intersect(rownames(expr), genes)
  print(paste("Num. of genes:", length(shared.genes)))
  expr <- expr[shared.genes, ]
  if (dim(expr)[1] < 5000) {
    print(paste("ERROR: not enough genes"))
    return - 1
  }
  
  # Transform the expression to rank
  expr <- apply(expr, 2, rank)
  
  # Run ssGSEA analysis for the ranked gene expression dataset
  if(packageVersion("GSVA") >= "1.36.0") {
    # GSVA >= 1.36.0 does not support `parallel.type` any more. 
    # Instead it automatically uses the backend registered by BiocParallel. 
    scores <- GSVA::gsva(expr, signatures, method = "ssgsea",
                         ssgsea.norm = FALSE,parallel.sz = parallel.sz)
  } else {
    scores <- GSVA::gsva(expr, signatures, method = "ssgsea",
                         ssgsea.norm = FALSE,parallel.sz = parallel.sz,parallel.type = parallel.type)
  }
  
  scores = scores - apply(scores,1,min)
  
  # Combine signatures for same cell types
  cell_types <- unlist(strsplit(rownames(scores), "%"))
  cell_types <- cell_types[seq(1, length(cell_types), 3)]
  agg <- aggregate(scores ~ cell_types, FUN = mean)
  rownames(agg) <- agg[, 1]
  scores <- agg[, -1]
  
  # Save raw scores
  if (!is.null(file.name)) {
    write.table(scores, file = file.name, sep = "\t",
                col.names = NA, quote = FALSE)
  }
  scores
}

genes <- xCell.data$genes
signatures <- xCell.data$signatures

row.names(tpm_denormalization) <- tpm_denormalization[,1]
tpm_denormalization <- tpm_denormalization[,-1]

rawscore <- rawEnrichmentAnalysis(tpm_denormalization, signatures, genes, 
                                  file.name = "xcell_rawscore_test.csv", 
                                  parallel.sz = 4, parallel.type = 'SOCK')


# To obatin immune microenvironment score
xCellAnalysis <- function(expr, signatures=NULL, genes=NULL, spill=NULL, rnaseq=TRUE, 
                          file.name = NULL, scale=TRUE,alpha = 0.5, save.raw = FALSE, 
                          parallel.sz = 4, parallel.type = 'SOCK',cell.types.use = NULL) {
  if (is.null(signatures))
    signatures = xCell.data$signatures
  if (is.null(genes))
    genes = xCell.data$genes
  if (is.null(spill)) {
    if (rnaseq==TRUE) {
      spill = xCell.data$spill
    } else {
      spill = xCell.data$spill.array
    }
  }
  
  # Caulcate average ssGSEA scores for cell types
  if (is.null(file.name) || save.raw==FALSE) {
    fn <- NULL
  } else {
    fn <- paste0(file.name,'_RAW.txt')
  }
  
  if (!is.null(cell.types.use)) {
    A = intersect(cell.types.use,rownames(spill$K))
    if (length(A)<length(cell.types.use)) {
      return ('ERROR - not all cell types listed are available')
    }
  }
  
  scores <- rawEnrichmentAnalysis(expr,signatures,genes,fn, parallel.sz = parallel.sz, 
                                  parallel.type = 'SOCK')
  
  # Transform scores from raw to percentages
  scores.transformed <- transformScores(scores, spill$fv, scale)
  
  # Adjust scores using the spill over compensation matrix
  if (is.null(file.name)) {
    fn <- NULL
  } else {
    fn <- file.name
  }
  
  if (is.null(cell.types.use)) {
    scores.adjusted <- spillOver(scores.transformed, spill$K, alpha,fn )
    scores.adjusted = microenvironmentScores(scores.adjusted)
  } else {
    scores.adjusted <- spillOver(scores.transformed[cell.types.use,], spill$K, alpha,fn )
  }
  return(scores.adjusted)
}

adj_df <- xCellAnalysis(tpm_denormalization)
adj_df <- data.frame(adj_df)


# Consensus clustering
set.seed(6)
GeneExp <- as.data.frame(rawscore)
GeneExp_m <- as.matrix(GeneExp)

result=ExecuteCC(clusterNum=3, GeneExp_m, reps=1000, maxK=5,clusterAlg="pam",distance="pearson",
                 title="tpm_consensus_3_cluster_reps_pearson", writeTable = TRUE)

clustered_group <- as.data.frame(t(result$group))
rownames(clustered_group) <- c("xCell")

xCell_clustered <- rbind(GeneExp, clustered_group) # gene expression data & clustered_group merged

xCell_clustered_t <- data.frame(t(xCell_clustered))