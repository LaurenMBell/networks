# ========================== CSF ============================================
"""
STEP 1: 
     - Raw data in Box, PROJECTS.Brain-gut_SaraCarloni.Metabolomics
     - Drop samples according to QC
     - Filter for CTRLs
     - Log transformation
     - Quantile normalization
     - Drop 0s for NaNs
"""

import pandas as pd
import numpy as np
import qnorm


def log_1(x):
    return np.log2(x + 1)


#TO KEEP
keep_samples = [
    "DSS_76_090224",
    "DSS_B896_090224",
    "DSS_T0.7_090224",
    "DSS_T0.8_090224",
    "DSS_73_090224",
    "935_CTR_LPS_CSF_070224",
    "936_CTR_LPS_CSF_070224_2",
    "937_CTR_LPS_CSF_070224",
    "941_CTR_LPS_CSF_070224",
    "942_CTR_LPS_CSF_070224",
    "943_CTR_LPS_CSF_070224",
    "3297_VECPAC_190224",
    "3679_VECPAC_190224",
    "3681_VECPAC_190224"]

df = pd.read_csv("CSF_data/CSF_raw.csv", header=None)

samples = df.iloc[0, 2:].str.strip().str.replace(r'\s+', '_', regex=True).values
metabolites = df.iloc[2:, 0].values
data = df.iloc[2:, 2:].values

data = pd.DataFrame(data, index=metabolites, columns=samples)
data = data.apply(pd.to_numeric, errors="coerce")

data = data[[col for col in data.columns if col in keep_samples]]

log_transformed = data.apply(log_1)

qnormed = qnorm.quantile_normalize(log_transformed, axis=1)

qnormed_nans = qnormed.copy()
qnormed_nans[data == 0] = np.nan

qnormed_nans.to_csv("CSF_preprocessed.csv", index=True)


#====================== CPX =================================

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

#filter cpx data to just controls
dss = pd.read_csv("CPX_data/DSS.csv")
vecpac = pd.read_csv("CPX_data/VECPAC.csv")
lps = pd.read_csv("CPX_data/LPS.csv") 

dss_map = pd.read_csv("CPX_data/DSS_group_map.csv", header=None)
vecpac_map = pd.read_csv("CPX_data/VECPAC_group_map.csv", header=None)
lps_map = pd.read_csv("CPX_data/LPS_group_map.csv", header=None)

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

mouse_map = pd.read_csv("CPX_data/mouse_map.csv")
mouse_map.columns = ["matt_id", "sara_id"] 
id_map = dict(zip(mouse_map["matt_id"].astype(str), mouse_map["sara_id"].astype(str)))
id_list = mouse_map["sara_id"].astype(str).tolist()

#change the column names of the cpx mice
cpx_dss = rename_cpx_columns(cpx_dss, mouse_map)
cpx_vecpac = rename_cpx_columns(cpx_vecpac, mouse_map)
cpx_lps = rename_cpx_columns(cpx_lps, mouse_map)


#=============== JOIN BY MOUSE ==========================
# Parse CSF sample names to extract mouse IDs for matching
csf_dss = qnormed_nans[[col for col in qnormed_nans.columns if 'DSS' in col]]
csf_vecpac = qnormed_nans[[col for col in qnormed_nans.columns if 'VECPAC' in col]]
csf_lps = qnormed_nans[[col for col in qnormed_nans.columns if 'LPS' in col]]

common_dss = sorted(list(set(cpx_dss.columns) & set(csf_dss.columns)))
common_vecpac = sorted(list(set(cpx_vecpac.columns) & set(csf_vecpac.columns)))
common_lps = sorted(list(set(cpx_lps.columns) & set(csf_lps.columns)))

# simplify to the mice they have in common
cpx_aligned_dss = cpx_dss[common_dss]
cpx_aligned_vecpac = cpx_vecpac[common_vecpac]
cpx_aligned_lps = cpx_lps[common_lps]
csf_aligned_dss = csf_dss[common_dss]
csf_aligned_vecpac = csf_vecpac[common_vecpac]
csf_aligned_lps = csf_lps[common_lps]

# merge each model by mouse
merged_dss = pd.concat([cpx_aligned_dss, csf_aligned_dss], axis=0)
merged_vecpac = pd.concat([cpx_aligned_vecpac, csf_aligned_vecpac], axis=0)
merged_lps = pd.concat([cpx_aligned_lps, csf_aligned_lps], axis=0)

print("MERGED DSS: ", merged_dss.columns)
merged_dss.to_csv("merged_dss.csv", index=False)
print("MERGED VECPAC: ", merged_vecpac.columns)
merged_vecpac.to_csv("merged_vecpac.csv", index=False)
print("MERGED LPS: ", merged_lps.columns)
merged_lps.to_csv("merged_lps.csv", index=False)


