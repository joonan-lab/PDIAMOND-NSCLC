#### Code title: run_singlecellAtlas
#### Description: Single-cell distribution on five subtypes using Single-cell atlas
#### Last log: 23.07.14

## Settings
rm(list = ls())

## Needed libraries
library(Seurat)
library(dplyr)
library(rtracklayer)
library(stringr)
library(ggplot2)
library(cowplot)
library(RColorBrewer)

## Load single cell atlas data
seurat.obj = readRDS("local_extended.rds")
seurat.obj.tp = subset(seurat.obj, subset = origin == "tumor_primary")
seurat.obj.tp <- SetIdent(seurat.obj.tp, value = "cell_type_tumor")
seurat.obj.tp.nat = subset(seurat.obj, subset = origin %in% c("tumor_primary", "normal_adjacent"))
seurat.obj.tp.nat <- SetIdent(seurat.obj.tp.nat, value = "cell_type_major")

## Load differential expression analysis result
rna.subtype = readRDS("DEG_DESeq2_bySubtype.Rdata") # 'run_DEG.R' script's output file is used for this input
rna.nat = readRDS("files/DEG_DESeq_byNAT.Rdata") # 'run_DEG.R' script's output file is used for this input

## Run (subtype and other subtypes)
rna.subtype.list = list()
for (cls in names(rna.nat)){
  rna.tmp = rna.subtype[[cls]]
  rna.tmp.sig = rna.tmp %>% filter(padj < 0.05 & log2FoldChange > 0) %>% filter(!is.na(gene_id))
  rna.tmp.sig$gene_id = sapply(strsplit(as.character(rna.tmp.sig$gene_id), "[.]"), "[", 1)
  rna.tmp.sig = rna.tmp.sig %>% pull("gene_id") %>% unique
  rna.subtype.list[[cls]] <- rna.tmp.sig
}

seurat.obj.tp <- AddModuleScore(seurat.obj.tp,
                                features = rna.subtype.list,
                                name="DEG_Subtype")

# More detail in neutrophil
seurat.type = seurat.obj.tp
meta = seurat.type@meta.data
meta2 = meta %>% mutate(cell_type_major_neutro = ifelse(grepl('NAN-', cell_type_neutro), 'NAN',
                                                        ifelse(grepl('TAN-', cell_type_neutro), 'TAN', as.character(cell_type_major))))

seurat.type@meta.data <- meta2
seurat.type <- SetIdent(seurat.type, value = "cell_type_major_neutro")

# UMAP draw
p.subtype = list()
for (feature in c("DEG_Subtype1", "DEG_Subtype2", "DEG_Subtype3", "DEG_Subtype4", "DEG_Subtype5")){
  p <- FeaturePlot(seurat.type, features = feature, label = TRUE, repel = TRUE, raster=TRUE) +
    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
  p.subtype[[feature]] <- p
}
plot_grid(p.subtype[[1]], p.subtype[[2]], p.subtype[[3]],  p.subtype[[4]], p.subtype[[5]], nrow= 1)


### Run (Tumor versus NAT on each subtype)
rna.nat.list = list()
for (cls in names(rna.nat)){
  rna.tmp = rna.nat[[cls]]
  rna.tmp.sig = rna.tmp %>% filter(padj < 0.05 & log2FoldChange > 0) %>% filter(!is.na(gene_id))
  rna.tmp.sig$gene_id = sapply(strsplit(as.character(rna.tmp.sig$gene_id), "[.]"), "[", 1)
  rna.tmp.sig = rna.tmp.sig %>% pull("gene_id") %>% unique
  rna.nat.list[[cls]] <- rna.tmp.sig
}

seurat.obj.tp.nat <- AddModuleScore(seurat.obj.tp.nat,
                                    features = rna.nat.list,
                                    name="DEG_NAT_Subtype")

# UMAP draw
p.nat = list()
for (feature in c("DEG_NAT_Subtype1", "DEG_NAT_Subtype2", "DEG_NAT_Subtype3", "DEG_NAT_Subtype4", "DEG_NAT_Subtype5")){
  p <- FeaturePlot(seurat.obj.tp.nat, features = feature, label = TRUE, repel = TRUE, raster=TRUE) +
    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
  p.nat[[feature]] <- p
}
plot_grid(p.nat[[1]], p.nat[[2]], p.nat[[3]], p.nat[[4]], p.nat[[5]], nrow = 1)

