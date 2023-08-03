#### Code title: run_facets_inferWGD
#### Description: WGD inference using FACETS
#### Last log: 23.07.24

## Settings
rm(list = ls())

## Needed libraries
library(tidyverse)
library(pctGCdata)
library(facets)

## Load SNP matrix (in csv file) for each sample called from cnv_facets.R (https://github.com/dariober/cnv_facets)
## The csv file contains 12 columns in total: Chromosome, Position, Ref, Alt, File1R, File1A, File1E, File1D, File2R, File2A, File2E, and File2D
fs = list.files('files/cnv_facets_results/', pattern = 'csv.gz')

# Call CNV with facets and merge output from all samples (https://github.com/mskcc/facets)
fit.list = list()
res = data.frame()
for (f in fs){
  print (f)
  sid = gsub('.csv.gz', '', strsplit(f, '_')[[1]][2])
  set.seed(1234)
  rcmat = readSnpMatrix(paste('files/cnv_facets_results/', f, sep=''))
  xx = preProcSample(rcmat, gbuild = "hg38")
  oo=procSample(xx, cval=150)
  fit.list[[sid]] <- emcncf(oo)
  res = rbind.data.frame(res, data.frame(sid = sid, purity=fit.list[[sid]]$purity, 
                                         ploidy=fit.list[[sid]]$ploidy, dipLogR=fit.list[[sid]]$dipLogR, emflags=ifelse(is.null(fit.list[[sid]]$emflags), 'None', fit.list[[sid]]$emflags)))
}

# Infer WGD (https://github.com/vanallenlab/facets)
wgd.res = data.frame()
for (s in names(res)){
  print (s)
  write.table(res[[s]]$cncf, 'in.txt', sep='\t', quote = F, row.names = F, col.names = T)
  system(command = 'python infer_wgd.py --cncf in.txt')
  
  d1 = read.delim('facets_wgd_boolean.txt', header = F)
  d2 = read.delim('facets_fraction_mcn_ge2.txt', header = F)
  
  wgd.res = rbind.data.frame(wgd.res, data.frame(sid=s, facet_wgd=d1$V1, facet_wgd_fraction=d2$V1))
  
  system(command = 'rm in.txt')
}

wgd.res = wgd.res %>% mutate(facet_wgd_fraction = facet_wgd_fraction * 100) 