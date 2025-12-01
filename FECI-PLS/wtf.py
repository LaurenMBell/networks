import pandas as pd
table1 = pd.read_csv("PLS-CPX_FDR_matt.csv")
table2 = pd.read_csv("PLS-CPX_FDR_reordered.csv")

table1.n1 = table1.n1.str.strip("-P")

table1_gene = set(table1['n1'])
table2_gene = set(table2['n1'])

print(len(table1_gene.difference(table2_gene))) 
#print(table1_gene.difference(table2_gene))

table1_met = set(table1['n2'])
table2_met = set(table2['n2'])

print(len(table1_met.difference(table2_met)))

pls_dss = pd.read_csv("metabolomics_data/pls_DSS_post_norm.csv")
pls_vecpac= pd.read_csv("metabolomics_data/pls_VECPAC_post_norm.csv")
pls_lps = pd.read_csv("metabolomics_data/pls_LPS_post_norm.csv")

cpx_dss = pd.read_csv("transcriptomics_data/DSS.csv")
cpx_vecpac = pd.read_csv("transcriptomics_data/VECPAC.csv")
cpx_lps = pd.read_csv("transcriptomics_data/LPS.csv")

og_genes = set(cpx_dss['ID'] + cpx_vecpac['ID'] + cpx_lps['ID'])
og_mets = set(pls_dss['ID'] + pls_vecpac['ID'] + pls_lps['ID'])

print(len(og_genes))





