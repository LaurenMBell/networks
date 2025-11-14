import pandas as pd
import numpy as np

"""
What you're going to do:
1) Filter to control samples.
2) Combine PLS and CPX data tables (join by mouse).
3) Compute pooled correlations (for p-value, correlation coefficient)
4) Compute grouped correlation (for correlation coefficient) (within VECPAC, LPS, & DSS)
5) Filter to edges where sign(pooled) = sign(VECPAC) = sign(LPS) = sign(DSS)
6) FDR correct correlation p-values
"""

def rename_pls_columns(df, id_list):
    new_cols = []
    
    #take the first token that matches the mouse_map, and the second if not
    for col in df.columns:
        tokens = col.split("_")
        second = tokens[1] if len(tokens) > 1 else col
        
        matched = False
        for id_str in id_list:
            if id_str in tokens:
                new_cols.append(id_str)
                matched = True
                break
        
        if not matched:
            new_cols.append(second)
    
    df.columns = new_cols
    return df

def rename_cpx_columns(df, id_map):
    df.columns = [id_map.get(str(col), col) for col in df.columns]

 #============== FILTER CPX and PLS DATA ====================

#filter cpx data to just controls
dss = pd.read_csv("transcriptomics_data/DSS.csv")
vecpac = pd.read_csv("transcriptomics_data/VECPAC.csv")
lps = pd.read_csv("transcriptomics_data/LPS.csv") 

dss_map = pd.read_csv("transcriptomics_data/DSS_group_map.csv", header=None)
vecpac_map = pd.read_csv("transcriptomics_data/VECPAC_group_map.csv", header=None)
lps_map = pd.read_csv("transcriptomics_data/LPS_group_map.csv", header=None)

#get just the control samples from each map
dss_map_controls = dss_map[dss_map[1] == "control"] 
vecpac_map_controls = vecpac_map[vecpac_map[1] == "control"]
lps_map_controls = lps_map[lps_map[1] == "control"]

#add ID col and new filtered data to get control samples for each model
dss_filtered = dss[["ID"] + dss_map_controls[0].astype(str).to_list()]
vecpac_filtered = vecpac[["ID"] + vecpac_map_controls[0].astype(str).to_list()]
lps_filtered = lps[["ID"] + lps_map_controls[0].astype(str).to_list()]

cpx_dss = dss_filtered.copy()
cpx_vecpac = vecpac_filtered.copy()
cpx_lps = lps_filtered.copy()

pls_dss = pd.read_csv("metabolomics_data/pls_DSS_post_norm.csv")
pls_vecpac = pd.read_csv("metabolomics_data/pls_VECPAC_post_norm.csv")
pls_lps = pd.read_csv("metabolomics_data/pls_LPS_post_norm.csv")

#============= CHANGE NAMES OF BOTH DATASETS ===========================

mouse_map = pd.read_csv("transcriptomics_data/mouse_map.csv")
mouse_map.columns = ["num_id", "sara_id"]   # adjust if your file already has headers
id_map = dict(zip(mouse_map["num_id"].astype(str), mouse_map["sara_id"].astype(str)))
id_list = mouse_map["sara_id"].astype(str).tolist()

#change the column names of the cpx mice
for df in [cpx_dss, cpx_vecpac, cpx_lps]:
    rename_cpx_columns(df, id_map)

#change the column names of the pls mice
for df in [pls_dss, pls_vecpac, pls_lps]:
    rename_pls_columns(df, id_list)

print(pls_dss.columns)
print(cpx_dss.columns)
print(pls_vecpac.columns)
print(cpx_vecpac.columns)
print(pls_lps.columns)
print(cpx_lps.columns)

#=============== JOIN BY MOUSE ==========================
common_dss = list(set(cpx_dss.columns) & set(pls_dss.columns))
common_vecpac = list(set(cpx_vecpac.columns) & set(pls_vecpac.columns))
common_lps = list(set(cpx_lps.columns) & set(pls_lps.columns))

cpx_aligned_dss = cpx_dss[common_dss]
cpx_aligned_vecpac = cpx_vecpac[common_vecpac]
cpx_aligned_lps = cpx_lps[common_lps]

pls_aligned_dss = pls_dss[common_dss]
pls_aligned_vecpac = pls_vecpac[common_vecpac]
pls_aligned_lps = pls_lps[common_lps]

merged_dss = pd.concat([cpx_aligned_dss, pls_aligned_dss], axis=0)
merged_vecpac = pd.concat([cpx_aligned_vecpac, pls_aligned_vecpac], axis=0)
merged_lps = pd.concat([cpx_aligned_lps, pls_aligned_lps], axis=0)

print(merged_dss.columns)
print(merged_vecpac.columns)
print(merged_lps.columns)

#================ POOLED CORRELATIONS ====================
