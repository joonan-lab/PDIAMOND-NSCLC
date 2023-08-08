#### Code title: run_MHCflurry
#### Description: Binding prediction of amino acid sequence and HLA allele
#### Last log: 23.08.04

## Needed libraries
from mhcflurry import Class1PresentationPredictor

## Load data
predictor = Class1PresentationPredictor.load()

## Load sample input file

# You can make by processing from maf file "Output_s230_20210827_Fun38_filtered_lite.Mutect2_maftools.maf" (maf -> vcf -> Ensembl VEP -> pVACtools)
# "Output_s230_20210827_Fun38_filtered_lite.Mutect2_maftools.maf": maf format file
# “RE-P001_20210827_Fun38_9mer.txt”: A file containing peptide fragments of a fixed length (9-12 amino acid sequence) on each line
input_sequence_file='RE-P001_20210827_Fun38_9mer.txt'

input_hla_allele='HLA-A02:06'

## Run MHCflurry

for line in open(fpath): #Read sequence file
    seqs+=[line.rstrip()]
    
result=predictor.predict(
    peptides=seqs,
    alleles=[allele],
    verbose=0)
cut_result=result.drop(['peptide_num','sample_name','best_allele'],axis=1)