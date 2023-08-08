#### Code title: run_NetMHCpan
#### Description: Binding prediction of amino acid sequence and HLA allele
#### Last log: 23.08.04

## Needed libraries
import os

## Load sample input file

# You can make by processing from maf file "Output_s230_20210827_Fun38_filtered_lite.Mutect2_maftools.maf" (maf -> vcf -> Ensembl VEP -> pVACtools)
# "Output_s230_20210827_Fun38_filtered_lite.Mutect2_maftools.maf": maf format file
# “RE-P001_20210827_Fun38_9mer.txt”: A file containing peptide fragments of a fixed length (9-12 amino acid sequence) on each line
input_sequence_file='RE-P001_20210827_Fun38_9mer.txt'

input_hla_allele='HLA-A02:06'

## Run NetMHCpan

NetMHCpan_path='/home/user/netMHCpan-4.1/netMHCpan' # your NetMHCpan path

os.system(f'{NetMHCpan_path} -p {input_sequence_file} -a {input_hla_allele} -BA > output_NetMHCpan.txt')
