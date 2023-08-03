# Proteogenomic Analysis of a Korean NSCLC Cohort

This study presents genetic proteomic data from 229 Korean patients with non-small cell lung cancer (NSCLC) to uncover the molecular characteristics of different subtypes, including those previously reported (LUAD/LSCC) and newly discovered ones. The paper serves as a valuable resource, providing an in-depth analysis of each subtype's molecular features.

We achieve this by conducting a comprehensive analysis of the subtypes' molecular biology, employing NMF subtyping, immune landscape assessment, Wgd phenotype investigation, survival analysis, and leveraging publicly available single-cell data.

## [1] Data resource

### Proteome

You can request proteomic raw data from K-BDS (Korea BioData Station, [https://kbds.re.kr](https://kbds.re.kr/)) with the accession ID KAP503860, KAP504787, KAP504788.

### Genome and Transcriptome

You can request genomic and transcriptomic raw data from Korean Nucleotide Archive (KoNA, [https://kobic.re.kr/kona](https://kobic.re.kr/kona)) with the accession ID, KAP210028.

## [2] Nonnegative Matrix Factorization (NMF) clustering

NMF clustering is a data analysis technique that identifies patterns and groups in data using Non-negative Matrix Factorization. For more information about this tool, see “[http://renozao.github.io/NMF/devel/index.html](http://renozao.github.io/NMF/devel/index.html)”

### Run NMF

We perform NMF clustering using expression data of global proteome, phospho proteome, and acetyl proteome. To get NMF results and view the heatmap, run the **"run_NMF.R"** script with your matrix for analysis. It will generate the NMF result file along with the heatmap image.

If you use your data, first follow the guidelines of the NMF tool package to analyze the optimal rank value. After analyzing with the optimal rank, you can visualize the pheatmap by inputting the data that corresponds to your samples’ information.

```r
# Check this parameter for optimal rank test in **“run_NMF.R”** script

n_run = 200 # If do optimal test, set 50.
mi = 5000 # If do optimal test, set 2000.
n_rank = 5 # Test of optimal rank for data should be first.
```

## [3] Feature-wise survival analysis

In order to understand the molecular basis of survival associations in certain subtypes, we performed a feature-wise survival analysis, identifying specific molecules associated with patient survival outcomes. This analysis helps uncover potential prognostic markers and therapeutic targets for personalized treatment strategies.

### Run survival analysis

Run **“run_Survival.R”** script with feature expression data.

We divided the samples into two groups, based on the expression levels of specific features, using the quantile() function to distinguish between high and low expression groups. We used cancer-specific overall survival length as the outcome variable to investigate how feature expression levels relate to patient survival durations.

## [4] Single-cell specific Subtype distribution analysis

We evaluated the tumor microenvironment of the five NSCLC subtypes by comparing subtype-specific genes with a diverse set of cell type-specific genes. This analysis was performed using an integrated single-cell RNA (scRNA) sequencing dataset from NSCLC patients ([https://luca.icbi.at](https://luca.icbi.at/)). 

### Run Differentially expressed (DE) analysis

Subtype-specific genes were identified by analyzing the transcriptome data using DESeq2. By using DESeq2 ([https://bioconductor.org/packages/release/bioc/html/DESeq2.html](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)), our analysis was conducted following these established guidelines to ensure reliable and robust results in detecting differentially expressed genes among the compared samples or subtypes.

To conduct the Tumor vs NAT (Normal Adjacent Tissue) analysis within each subtype and the Tumor of Subtype vs Tumor of Other Subtypes analysis, execute the “**run_DEG.R”** script. Keep in mind that the DESeq2 analysis design will differ for each analysis as outlined below.

```r
# 1. Tumor vs NAT analysis within each subtype:
design = ~ Group + rna.batch

# 2. Tumor of Subtype vs Tumor of Other Subtypes analysis:
design = ~ Group + DX + rna.batch
```

By running the “**run_DEG.R”** script with these specific designs, you can perform the differential gene expression analysis as required. The obtained results from the DESeq2 analysis were further used for single-cell analysis, applying the filtering criteria of 'padj (Adjusted p value) < 0.05' and 'log2FoldChange > 0'. This filtering process helps to focus on statistically significant genes with a significant fold change, which are likely to be biologically relevant and informative for the subsequent single-cell analysis.

### Run single-cell analysis

In our analysis, we utilized the file named "local_extended.rds", which is available for download from the open-source repository at **[https://luca.icbi.at](https://luca.icbi.at/)**. This resource is likely to provide additional data and extended information that complements our research and enhances the results. 

Run “**run_singlecellAtlas.R”** script with above DE analysis results. In our analysis focused on tumors in subtypes, we excluded the 'normal_adjacent' origin from the distribution analysis. Indeed, users have the flexibility to adjust the subset objects according to their specific research needs and purposes.

```r
# Load single cell atlas data
seurat.obj = readRDS("files/local_extended.rds")
seurat.obj.tp = subset(seurat.obj, subset = origin == "tumor_primary")
seurat.obj.tp <- SetIdent(seurat.obj.tp, value = "cell_type_tumor")
seurat.obj.tp.nat = subset(seurat.obj, subset = origin %in% c("tumor_primary", "normal_adjacent"))
seurat.obj.tp.nat <- SetIdent(seurat.obj.tp.nat, value = "cell_type_major")
```

After obtaining the results, we proceeded to visualize the UMAP (Uniform Manifold Approximation and Projection) to represent the subtype-specific cell distributions. 

## [5] Whole-genome doubling (WGD) analysis

### Inferring whole-genome doubling (WGD) status

To determine the WGD status in each sample, we employed a multi-step process. First, we generated SNP matrices using the cnv_facets tool ([https://github.com/dariober/cnv_facets](https://github.com/dariober/cnv_facets)). The input for this process consisted of preprocessed paired tumor-normal BAM files and a VCF file containing common and germline polymorphic sites, which were obtained from [https://www.ncbi.nlm.nih.gov/variation/docs/human_variation_vcf/](https://www.ncbi.nlm.nih.gov/variation/docs/human_variation_vcf/). Next, we utilized FACETS ([https://github.com/mskcc/facets](https://github.com/mskcc/facets)) to call copy number variations (CNVs). Subsequently, we inferred the WGD status ([https://github.com/vanallenlab/facets](https://github.com/vanallenlab/facets)) by adhering to established guidelines.

To replicate this analysis, first generate the SNP matrices using preprocessed paired tumor-normal BAM files, then execute the **"run_facets_inferWGD.R"** script.

### Copy number signature analysis

For the copy number signature analysis, we employed the COSMIC signature database v3 ([https://cancer.sanger.ac.uk/signatures/](https://cancer.sanger.ac.uk/signatures/)) and the R package Sigminer ([https://github.com/ShixiangWang/sigminer](https://github.com/ShixiangWang/sigminer)), following the instructions provided at [https://shixiangwang.github.io/sigminer-doc/cnsig.html](https://shixiangwang.github.io/sigminer-doc/cnsig.html).

As input data, we utilized the CNV data obtained from the previous analysis. We merged this data across samples and modified it to have six columns with the following column names: Chromosome, Start.bp, End.bp, modal_cn, minor_cn, and sample. To determine the number of signature groups (or factorization rank), we employed non-negative matrix factorization (NMF) with a tumor-by-component matrix and performed 50 runs while checking 2 to 12 ranks. Based on the cophenetic score plot, we chose rank seven for copy number variants. Each signature was identified using the COSMIC signature "CNS_TCGA" with the highest cosine similarity. We then assigned samples to one of the signatures based on the consensus matrix through hierarchical clustering.

To conduct copy number signature analysis with the same settings, you can execute the "**run_CNsignature.R**" script using your CNV input. Additionally, you have the flexibility to modify parameters and the reference signature database according to your specific requirements.

## [6] Tumor immune microenvironment analysis

Tumor immune microenvironment (TIME) plays a critical role in understanding the molecular characteristics of various cancer types. To gain insights into these molecular features, we conducted TIME profiling and identified genes that influence the patient's TIME.

### Cell composition inference and clustering of the tumor immune microenvironment

In our analysis, we employed two methods to cluster the tumor immune microenvironment. One method is xCell ([https://github.com/dviraran/xCell](https://github.com/dviraran/xCell)), which is an analysis based on cell type enrichment, and the other method involves performing GSVA analysis ([https://bioconductor.org/packages/release/bioc/html/GSVA.html](https://bioconductor.org/packages/release/bioc/html/GSVA.html)) based on immune-related pathways described [Hu et al](https://www.frontiersin.org/journals/oncology/articles/10.3389/fonc.2020.01189/full). Based on the cell composition socres obtained from xCell, we utilized the "CancerSubtypes" package in R ([https://bioconductor.org/packages/release/bioc/html/CancerSubtypes.html](https://bioconductor.org/packages/release/bioc/html/CancerSubtypes.html)) to perform clustering into hot-tumor-enriched (HTE), cold-tumor-enriched (CTE), and NAT-enriched subtypes.

```python
# Parameters for Execute Consensus Clustering
clusterNum=3
reps=1000
clusterAlg="pam"
distance="pearson"
```

### Identification of putative regulators associated with immune landscape

To identify key genes influencing the immune landscape, we used transcriptome data and proteomics data. We inferred protein activity based on the transcriptome data using the "viper" package in R ([https://www.bioconductor.org/packages/release/bioc/html/viper.html](https://www.bioconductor.org/packages/release/bioc/html/viper.html)). To utilize this tool, we downloaded the regulon network from the "ARACNe" package in R ([http://bioconductor.org/packages/release/data/experiment/html/aracne.networks.html](http://bioconductor.org/packages/release/data/experiment/html/aracne.networks.html)). Next, we examined the correlation between the levels of RNA expression, protein expression, and protein activity in relation to the enrichment score of immune cells or the status of immune clusters and immunomodulators using the "stats" package in R.

Furthermore, to analyze whether specific driver mutations are involved in the relationship between immunomodulators and immune cells or immune clusters, we utilized driver mutation data from Oncovar ([Wang et al](https://academic.oup.com/nar/article/49/D1/D1289/5976976?login=true)), DriverDBv3 ([Liu et al](https://academic.oup.com/nar/article/48/D1/D863/5614573?login=true)), Intogen ([Gundem et al](https://www.nature.com/articles/nmeth0210-92)), and the mutation catalogue from [Martínez-Jiménez et al](https://www.nature.com/articles/s41568-020-0290-x).

## [7] Neo-antigen prediction

### Prediction of HLA binding neoepitopes

To generate a peptide fragment from a mutation locus, we employed Ensembl VEP (Variant Effect Predictor) to obtain detailed information for each mutation in the VCF file derived from Whole Exome Sequencing (WES) ([https://asia.ensembl.org/info/docs/tools/vep/index.html](https://asia.ensembl.org/info/docs/tools/vep/index.html)). 

```powershell
# Parameters for vep run
--no_stats -vcf --force --symbol --terms SO --tsl --hgvs --offline --pick --transcript_version --no_check_variants_order --buffer_size 20000

# Plugin
--plugin Downstream --plugin Frameshift --plugin Wildtype

# Reference genome
--fasta Homo_sapiens.GRCh38.dna.primary_assembly.fa
```

Afterwards, we extended the flanking amino acid sequences around the mutations using pVACseq, which is a module of pVACtools ([https://pvactools.readthedocs.io/en/latest/](https://pvactools.readthedocs.io/en/latest/)). These amino acid sequences were segmented into lengths of 9-12 amino acids to be predict their binding to HLAs. This was achieved by ustilizing NetMHCpan ([https://services.healthtech.dtu.dk/services/NetMHCpan-4.1/](https://services.healthtech.dtu.dk/services/NetMHCpan-4.1/)) and MHCflurry ([https://github.com/openvax/mhcflurry](https://github.com/openvax/mhcflurry)) against patient-specific HLA alleles. In our consideration of HLA binding neoepitopes, we defined them as 'WB' (Weak Binder) or 'SB' (Strong Binder) based on the predictions from NetMHCpan. Additionally, we identified those with a binding affinity under 500nM in the predictions from MHCflurry.