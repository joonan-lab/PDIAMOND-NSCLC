#### Code title: run_DEG
#### Description: Differential expression analysis for all variables using this study
#### Last log: 23.07.14

## Settings
rm(list = ls())

## Needed libraries
library(tximport)
library(tidyverse)
library(DESeq2)

### Run DEG analysis (subtype vs other subtypes)
## Load data (You can download in supplementary file)
s = readxl::read_xlsx('TableS1', sheet = 2) %>% mutate(DX=as.character(DX))
s$sid = s$Sample.ID

fs = list.files('files/nf_star_salmon/')
samples <- data.frame(sid=gsub('RE-', 'RE-P', fs), 
                      file=paste('files/nf_star_salmon/', fs, '/quant.sf', sep='')) %>% 
  mutate(tissue=ifelse(grepl('-T',fs), 'T', 'N'), tissue=factor(tissue)) %>% filter(tissue=='T')
tx2gene = read_delim('files/salmon_tx2gene.tsv', col_names = F) %>% select(transcript_id = 1, gene_id = 2, gene_name = 3)
gm = rtracklayer::import('files/gencode.v32.annotation.gtf.gz')
gm1 = gm %>% as.data.frame() %>% filter(type=='transcript') %>% filter(transcript_type=='protein_coding') %>% 
  select(gene_name, gene_id) %>% unique()

## Add RNA STAR stats
d.log = readRDS('files/data_merged_nf_STAR_stats.20211113.Rdata')
d.log = d.log %>% filter(stat %in% c('Uniquely mapped', 'Uniquely mapped reads %')) %>% 
  pivot_wider(names_from = stat, values_from = value) %>% select(sid=1, rna.batch=2, 'unique_map_total'=3, 'unique_map_perc'=4) 
s1 = merge(s, d.log, by='sid', all.x=T)
s1 = s1 %>% mutate(rna.batch = as.factor(rna.batch))

## Subtype1 ~ Subtype5 (DEG analysis)
DDS = list()
deg_res = list()
deseq_counts = list()
for(to_test in 1:5){
  tag = paste("Subtype",to_test,sep="")
  print(tag)
  
  samples2 = merge(samples, s1 %>% select(sid, Pathologic.N, Subtype, rna.batch, DX, unique_map_perc, unique_map_total), by='sid') %>% 
    mutate(DX=factor(DX, levels = c('AD', 'SC', 'MA', 'NC', 'Others')), 
           rna.batch = gsub('-','',rna.batch),
           rna.batch = as.factor(rna.batch),
           Group = ifelse(Subtype==to_test, 'Y', 'N'), 
           Group = factor(Group, levels=c('N', 'Y')))
  
  ## Create a DESeq2 object
  txi = tximport(samples$file, type="salmon", tx2gene=tx2gene, txIn = T, txOut = F)
  dds = DESeqDataSetFromTximport(txi, colData = samples2, design = ~  Group + DX + rna.batch)
  dds = estimateSizeFactors(dds)
  
  DDS[[tag]] = dds
  
  # DEG analysis
  keep = rowMeans(counts(dds)) >= 50
  table(keep)
  dds = dds[keep,]
  
  dds = DESeq(dds, parallel = T) # do parallel
  resultsNames(dds) # lists the coefficients
  res = results(dds, name = "Group_Y_vs_N")
  
  res1 = data.frame(gene_id = rownames(res), baseMean=res$baseMean, log2FoldChange=res$log2FoldChange, lfcSE=res$lfcSE, stat=res$stat, pvalue=res$pvalue, padj=res$padj) %>% merge(., gm1, by='gene_id')
  
  deg_res[[tag]] = res1
}


### Analysis for NAT versus Tumor in each subtypes
s = readxl::read_xlsx('TableS1', sheet = 2) %>% mutate(DX=as.character(DX))
s$sid = s$Sample.ID

fs = list.files('files/nf_star_salmon/')
samples <- data.frame(sid=gsub('RE-', 'RE-P', fs), file=paste('files/nf_star_salmon/', fs, '/quant.sf', sep='')) %>% 
  mutate(tissue=ifelse(grepl('-T',fs), 'T', 'N'), tissue=factor(tissue))
tx2gene = read_delim('files/salmon_tx2gene.tsv', col_names = F) %>% 
  select(transcript_id = 1, gene_id = 2, gene_name = 3)
gm = rtracklayer::import('files/gencode.v32.annotation.gtf.gz')
gm1 = gm %>% as.data.frame() %>% filter(type=='transcript') %>% filter(transcript_type=='protein_coding') %>% 
  select(gene_name, gene_id) %>% unique()

# Sample annotation (Add NAT)
s.new = data.frame(sid = samples$sid)
s.new = s.new %>% merge(., s %>% dplyr::select(sid, DX, Subtype), by = 'sid', all.x=T)

## Add RNA STAR stats
d.log = readRDS('files/data_merged_nf_STAR_stats.20211113.Rdata')
d.log = d.log %>% filter(stat %in% c('Uniquely mapped', 'Uniquely mapped reads %')) %>% 
  pivot_wider(names_from = stat, values_from = value) %>% select(sid=1, rna.batch=2, 'unique_map_total'=3, 'unique_map_perc'=4) 
s1 = merge(s.new, d.log, by='sid', all.x=T)
s1 = s1 %>% mutate(rna.batch = as.factor(rna.batch))
s1 = s1 %>% mutate(Subtype = ifelse(is.na(Subtype), 'NAT', Subtype))
s1 = s1 %>% mutate(DX = ifelse(is.na(DX), 'NAT', DX))

## Subtype1 ~ Subtype5 (DEG analysis)
DDS = list()
deg_res = list()
deseq_counts = list()
for(to_test in 1:5){
  tag = paste("Subtype",to_test,sep="")
  print(tag)
  
  samples2 = merge(samples, s1 %>% select(sid, Subtype, rna.batch, DX, unique_map_perc, unique_map_total), by='sid') %>% 
    filter(Subtype %in% c(to_test, 'NAT')) %>%
    mutate(DX=factor(DX, levels = c('AD', 'SC', 'MA', 'NC', 'Others','NAT')), 
           rna.batch = gsub('-','',rna.batch),
           rna.batch = as.factor(rna.batch),
           Group = ifelse(Subtype==to_test, 'Y', 'N'), 
           Group = factor(Group, levels=c('N', 'Y')))
  
  ## Create a DESeq2 object
  txi = tximport(samples2$file, type="salmon", tx2gene=tx2gene, txIn = T, txOut = F)
  dds = DESeqDataSetFromTximport(txi, colData = samples2, design = ~  Group + rna.batch)
  dds = estimateSizeFactors(dds)
  
  DDS[[tag]] = dds
  
  # DEG analysis
  keep = rowMeans(counts(dds)) >= 50
  table(keep)
  dds = dds[keep,]
  
  dds = DESeq(dds, parallel = T) # do parallel
  resultsNames(dds) # lists the coefficients
  res = results(dds, name = "Group_Y_vs_N")
  
  res1 = data.frame(gene_id = rownames(res), baseMean=res$baseMean, 
                    log2FoldChange=res$log2FoldChange, lfcSE=res$lfcSE, 
                    stat=res$stat, pvalue=res$pvalue, padj=res$padj) %>% 
    merge(., gm1, by='gene_id')
  
  deg_res[[tag]] = res1
}

