#### Code title: run_Ensembl_VEP
#### Description: Annotating additional information about mutations
#### Last log: 23.08.04

## Needed libraries
import os

## Load sample input file

# You can make by processing from maf file "Output_s230_20210827_Fun38_filtered_lite.Mutect2_maftools.maf" (maf -> vcf)
# "Output_s230_20210827_Fun38_filtered_lite.Mutect2_maftools.maf": maf format file
# “Output_RE-P001_20210827_Fun38_filtered_lite.Mutect2_maftools.vcf”: vcf format file
input_file='Output_RE-P001_20210827_Fun38_filtered_lite.Mutect2_maftools.vcf'

# You can downloaded from "https://ftp.ensembl.org/pub/release-109/fasta/homo_sapiens/dna/"
reference_fasta='Homo_sapiens.GRCh38.dna.primary_assembly.fa'

## Run Ensembl VEP

vep_path='/home/user/ensembl-vep/vep' # your Ensembl VEP path
cache_path='/home/user/.vep/' # your cache directory path

os.system(f'{vep_path} --no_stats -vcf --force --symbol --terms SO --tsl --hgvs --fasta {reference_seq_fasta} --offline --pick --transcript_version --cache --plugin Downstream --plugin Frameshift --plugin Wildtype --no_check_variants_order --buffer_size 20000 --dir_cache {cache_path} --assembly GRCh38 -i {input_file} -o output_Ensembl_vep.vep.vcf')