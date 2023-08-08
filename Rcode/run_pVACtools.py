#### Code title: run_pVACtools
#### Description: Extending the front and back of an amino acid sequence
#### Last log: 23.08.04

## Needed libraries
import os

## Load sample input file

# You can make by processing from maf file "Output_s230_20210827_Fun38_filtered_lite.Mutect2_maftools.maf" (maf -> vcf -> Ensembl VEP)
# "Output_s230_20210827_Fun38_filtered_lite.Mutect2_maftools.maf": maf format file
# “Output_RE-P001_20210827_Fun38_filtered_lite.Mutect2_maftools.vep.vcf”: Ensembl VEP annotated vcf format file
input_file='Output_RE-P001_20210827_Fun38_filtered_lite.Mutect2_maftools.vep.vcf'

extend_length=8 # Length you want to extend

## Run pVACtools

os.system(f'pvacseq generate_protein_fasta {input_file} {extend_length} output_pVACtools.txt')