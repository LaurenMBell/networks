import pandas as pd
import numpy as np
from scipy import stats

def rename_feci_columns(df, id_list):
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

def rename_cpx_columns(df, mouse_map):
    num_to_sara = dict(zip(mouse_map["matt_id"].astype(str),
                           mouse_map["sara_id"].astype(str)))

    new_cols = []
    for col in df.columns.astype(str):
        if col == "ID":
            new_cols.append("ID")
            continue
        new_cols.append(num_to_sara[col])
    df.columns = new_cols
    return df

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
vecpac_map_controls = vecpac_map[vecpac_map[1] == "treatment"]
lps_map_controls = lps_map[lps_map[1] == "control"]

#add ID col and new filtered data to get control samples for each model
dss_filtered = dss[["ID"] + dss_map_controls[0].astype(str).to_list()]
vecpac_filtered = vecpac[["ID"] + vecpac_map_controls[0].astype(str).to_list()]
lps_filtered = lps[["ID"] + lps_map_controls[0].astype(str).to_list()]

cpx_dss = dss_filtered.copy()
cpx_vecpac = vecpac_filtered.copy()
cpx_lps = lps_filtered.copy()

feci_dss = pd.read_csv("feci_data/feci_DSS_post_norm.csv")
feci_vecpac = pd.read_csv("feci_data/feci_VECPAC_post_norm.csv")
feci_lps = pd.read_csv("feci_data/feci_LPS_post_norm.csv")

#============= CHANGE NAMES OF BOTH DATASETS ===========================

mouse_map = pd.read_csv("transcriptomics_data/mouse_map.csv")
mouse_map.columns = ["matt_id", "sara_id"] 
id_map = dict(zip(mouse_map["matt_id"].astype(str), mouse_map["sara_id"].astype(str)))
id_list = mouse_map["sara_id"].astype(str).tolist()

#change the column names of the cpx mice
cpx_dss = rename_cpx_columns(cpx_dss, mouse_map)
cpx_vecpac = rename_cpx_columns(cpx_vecpac, mouse_map)
cpx_lps = rename_cpx_columns(cpx_lps, mouse_map)

#change the column names of the feci mice
feci_dss = rename_feci_columns(feci_dss, id_list)
feci_vecpac = rename_feci_columns(feci_vecpac, id_list)
feci_lps = rename_feci_columns(feci_lps, id_list)

print("FECI - DSS:", feci_dss.columns)
print("CPX - DSS:",cpx_dss.columns)
print("FECI - VECPAC:",feci_vecpac.columns)
print("CPX - VECPAC:",cpx_vecpac.columns)
print("FECI - LPS:",feci_lps.columns)
print("CPX - LPS:", cpx_lps.columns)

#get rid of the control row
feci_dss = feci_dss.drop(0)
feci_vecpac = feci_vecpac.drop(0)
feci_lps = feci_lps.drop(0)

#=============== JOIN BY MOUSE ==========================
common_dss = sorted(list(set(cpx_dss.columns) & set(feci_dss.columns)))
common_vecpac = sorted(list(set(cpx_vecpac.columns) & set(feci_vecpac.columns)))
common_lps = sorted(list(set(cpx_lps.columns) & set(feci_lps.columns)))


# simplify to the mice they have in common
cpx_aligned_dss = cpx_dss[common_dss]
cpx_aligned_vecpac = cpx_vecpac[common_vecpac]
cpx_aligned_lps = cpx_lps[common_lps]
feci_aligned_dss = feci_dss[common_dss]
feci_aligned_vecpac = feci_vecpac[common_vecpac]
feci_aligned_lps = feci_lps[common_lps]

# merge each model by mouse
merged_dss = pd.concat([cpx_aligned_dss, feci_aligned_dss], axis=0)
merged_vecpac = pd.concat([cpx_aligned_vecpac, feci_aligned_vecpac], axis=0)
merged_lps = pd.concat([cpx_aligned_lps, feci_aligned_lps], axis=0)

print("MERGED DSS: ", merged_dss.columns)
merged_dss.to_csv("merged_dss.csv", index=False)
print("MERGED VECPAC: ", merged_vecpac.columns)
merged_vecpac.to_csv("merged_vecpac.csv", index=False)
print("MERGED LPS: ", merged_lps.columns)
merged_lps.to_csv("merged_lps.csv", index=False)



