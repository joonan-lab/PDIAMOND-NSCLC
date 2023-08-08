#### Code title: Confirm_regression_result_on_immunomodulator
#### Description: Investigating the impact of immunomodulators on immune cells or the status of immune clusters
#### Last log: 23.08.06

## Needed libraries
import sys, os
import pandas as pd
from tqdm import tqdm

# Load data
im_gene = pd.read_csv("Immunomodulator_gene_list.csv", index_col=0)
# immunomodulator gene list was download from Thorsson et al. "The immune landscape of cancer."
# "Immunomodulator_gene_list.csv"; Row: Immunomodulator category x Column: gene symbol.

r_df = pd.read_csv("RNA_expression_regression_result_pval05.csv", index_col=0)
pe_df = pd.read_csv("Protein_expression_regression_result_pval05", index_col=0)
pa_df = pd.read_csv("Protein_activity_regression_result_pval05", index_col=0)
# The above three datasets are the results of regression obtained through "run_immunomodulator_regression.R," using a p-value cutoff of 0.05.
# Three dataset formats; Row: gene symbol x Column: 11 major immune or stromal cells

rna_all_df = pd.read_csv("RNA_expression_regression_result.csv", index_col=0)
pe_all_df = pd.read_csv("Protein_expression_regression_result", index_col=0)
pa_all_df = pd.read_csv("Protein_activity_expression_regression_result.csv", index_col=0)
# The above three datasets are the regression results obtained through "run_immunomodulator_regression.R," encompassing the complete results.
# Three dataset formats; Row: gene symbol x Column: 11 major immune or stromal cells

# multiomics_data
g_r = pd.read_csv("gillette_RNA_expression_regression_result.csv", index_col=0)
g_pe = pd.read_csv("gillette_Protein_expression_regression_result.csv", index_col=0)
g_pa = pd.read_csv("gillette_Protein_activity_regression_result.csv", index_col=0)
s_r = pd.read_csv("satpathy_RNA_expression_regression_result.csv", index_col=0)
s_pe = pd.read_csv("satpathy_Protein_expression_regression_result.csv", index_col=0)
s_pa = pd.read_csv("satpathy_Protein_activity_expression_regression_result.csv", index_col=0)
# The above six datasets are the results obtained by performing regression on the expression values from "Satpathy et al" and "Gillette et al" using the same method as in "run_immunomodulator_regression.R."
# Six dataset formats; Row: gene symbol x Column: 11 major immune or stromal cells


# confirm regression result
r_df = r_df[r_df["xcell_subtype"]!=0]
im_gene = im_gene.set_index("gene")
im_gene["Data_type"] = "RNA"
im_gene = im_gene[["Data_type", "category"]]

r_f = im_gene.merge(r_df, left_index=True, right_index=True)
r_f = r_f.reset_index(drop=False)
r_f.columns = ["Gene_name", "Data_type", "Category", "CD8_T_cells", 
               "CD4_T_cells", "Tregs", "B_cells", "NK_cells", 
               "Neutrophils", "DC", "Monocytes", "Macrophages", 
               "Epithelial_cells", "Features"]
merge = pd.merge(xcell_cli, rna.T, left_index=True, right_index=True, how="left")


pe_df = pe_df[pe_df["xcell_subtype"]!=0]
im_gene["Data_type"] = "Protein_exp"

pe_f = im_gene.merge(pe_df, left_index=True, right_index=True)
pe_f = pe_f.reset_index(drop=False)
pe_f.columns = ["Gene_name", "Data_type", "Category", "CD8_T_cells", 
               "CD4_T_cells", "Tregs", "B_cells", "NK_cells", 
               "Neutrophils", "DC", "Monocytes", "Macrophages", 
               "Epithelial_cells", "Features"]


pa_df = pa_df[pa_df["xcell_subtype"]!=0]
im_gene["Data_type"] = "Protein_act"

pa_f = im_gene.merge(pa_df, left_index=True, right_index=True)
pa_f = pa_f.reset_index(drop=False)
pa_f.columns = ["Gene_name", "Data_type", "Category", "CD8_T_cells", 
               "CD4_T_cells", "Tregs", "B_cells", "NK_cells", 
               "Neutrophils", "DC", "Monocytes", "Macrophages", 
               "Epithelial_cells", "Features"]

m_df = pd.concat([r_f, pe_f, pa_f])

r_g = list(m_df[m_df["Data_type"]=="RNA"]["Gene_name"])
pe_g = list(m_df[m_df["Data_type"]=="Protein_exp"]["Gene_name"])
pa_g = list(m_df[m_df["Data_type"]=="Protein_act"]["Gene_name"])

a_l = list(set(r_g) & set(pe_g) & set(pa_g))
t_l = list(set(r_g) & set(pe_g))
t_l2 = list(set(pa_g) & set(pe_g))

a_l.extend(t_l)
a_l.extend(t_l2)
im_g = list(set(a_l))

def filtering_df(data, im_gene):
    tmp_df = data[data["genes"].isin(im_gene)]
    tmp_df = tmp_df.set_index("genes")
    
    return tmp_df

r_f = filtering_df(rna_all_df, im_g)
pe_f = filtering_df(pe_all_df, im_g)
pa_f = filtering_df(pa_all_df, im_g)

def filtering_cptac_df(cptac):
    cptac.index = cptac.index.map(lambda x: x.split("_")[1])
    
    return cptac

g_rf = filtering_cptac_df(g_r)
g_pef = filtering_cptac_df(g_pe)
g_paf = filtering_cptac_df(g_pa)

s_rf = filtering_cptac_df(s_r)
s_pef = filtering_cptac_df(s_pe)
s_paf = filtering_cptac_df(s_pa)

s_m = pd.concat([r_f, pe_f, pa_f, s_rf, s_pef, s_paf, g_rf, g_pef, g_paf], axis=1)

im_gene = pd.read_csv(Immunomodulator_gene_list)
im_gene = im_gene.set_index("gene")
im_gene.head()

m_df = pd.merge(s_m, im_gene, left_index=True, right_index=True, how="left")

# filtering
m_df = m_df.loc[['SLAMF7', 'CD40', 'VEGFA', 'ENTPD1', 'GZMA',
                'ICOS', 'SELP', 'BTLA', 'CCL5']]

m_df = m_df.sort_values(by=["category"])
m_df = m_df.reset_index(drop=False)
m_df.columns = ["Gene_ID", "PDIA_RNA", "PDIA_RNA_p", "PDIA_protein_exp", 
                "PDIA_protein_exp_p", "PDIA_protein_act", "PDIA_protein_act_p", 
                "Satpathy_RNA", "Satpathy_RNA_p", "Satpathy_protein_exp", 
                "Satpathy_protein_exp_p", "Satpathy_protein_act", "Satpathy_protein_act_p",
                "Gillette_RNA", "Gillette_RNA_p", "Gillette_protein_exp", "Gillette_protein_exp_p",
                "Gillette_protein_act", "Gillette_protein_act_p", "Category"]

m_df = m_df[["Gene_ID", "Category", "PDIA_RNA", "PDIA_protein_exp", "PDIA_protein_act", 
             "Satpathy_RNA", "Satpathy_protein_exp", "Satpathy_protein_act", "Gillette_RNA", 
             "Gillette_protein_exp", "Gillette_protein_act", "PDIA_RNA_p", "PDIA_protein_exp_p", 
             "PDIA_protein_act_p", "Satpathy_RNA_p", "Satpathy_protein_exp_p", "Satpathy_protein_act_p", 
             "Gillette_RNA_p", "Gillette_protein_exp_p", "Gillette_protein_act_p"]]