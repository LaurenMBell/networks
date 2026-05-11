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


def csf_names(sample_name, sara_ids):
    tokens = str(sample_name).replace("_", " ").split()
    for token in tokens:
        if token in sara_ids:
            return token

        #the one B sample for some reason
        if token.startswith("B") and token[1:] in sara_ids:
            return token[1:]

    return sample_name


def rename_csf_columns(df, sara_ids):
    df = df.copy()
    df.columns = [csf_names(col, sara_ids) for col in df.columns]
    return df


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

samples = df.iloc[0, 2:].str.strip().replace(r'\s+', '_', regex=True).values
metabolites = df.iloc[2:, 0].values
data = df.iloc[2:, 2:].values

data = pd.DataFrame(data, index=metabolites, columns=samples)
data = data.apply(pd.to_numeric, errors="coerce")

data = data[[col for col in data.columns if col in keep_samples]]

log_transformed = data.apply(log_1)

qnormed = qnorm.quantile_normalize(log_transformed, axis=1)

qnormed_nans = qnormed.copy()
qnormed_nans[data == 0] = np.nan
csf_preprocessed = pd.DataFrame(qnormed_nans)

csf_preprocessed.to_csv("CSF_preprocessed.csv", index=True)


#====================== STR =================================
def rename_str_columns(df, mouse_map):
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

# filter cpx data to just controls
dss = pd.read_csv("STR_data/DSS.csv")
vecpac = pd.read_csv("STR_data/VECPAC.csv")
lps = pd.read_csv("STR_data/LPS.csv")

dss_map = pd.read_csv("STR_data/DSS_group_map.csv", header=None)
vecpac_map = pd.read_csv("STR_data/VECPAC_group_map.csv", header=None)
lps_map = pd.read_csv("STR_data/LPS_group_map.csv", header=None)

# get just the control samples from each map
dss_map_controls = dss_map[dss_map[1] == "control"]
vecpac_map_controls = vecpac_map[vecpac_map[1] == "treatment"]
lps_map_controls = lps_map[lps_map[1] == "control"]

# add ID col and new filtered data to get control samples for each model
dss_filtered = dss[["ID"] + dss_map_controls[0].astype(str).to_list()]
vecpac_filtered = vecpac[["ID"] + vecpac_map_controls[0].astype(str).to_list()]
lps_filtered = lps[["ID"] + lps_map_controls[0].astype(str).to_list()]

str_dss = dss_filtered.copy()
str_vecpac = vecpac_filtered.copy()
str_lps = lps_filtered.copy()

mouse_map = pd.read_csv("STR_data/mouse_map.csv")
mouse_map.columns = ["matt_id", "sara_id"]
id_map = dict(zip(mouse_map["matt_id"].astype(str), mouse_map["sara_id"].astype(str)))
id_list = mouse_map["sara_id"].astype(str).tolist()

# change the column names of the cpx mice
str_dss = rename_str_columns(str_dss, mouse_map)
str_vecpac = rename_str_columns(str_vecpac, mouse_map)
str_lps = rename_str_columns(str_lps, mouse_map)

# rename CSF columns to Sara IDs for joining
sara_ids = set(mouse_map["sara_id"].astype(str))

#=============== JOIN BY MOUSE ==========================
# parse CSF sample names to get mouse IDs for matching
csf_dss = csf_preprocessed[[col for col in csf_preprocessed.columns if 'DSS' in col]]
csf_vecpac = csf_preprocessed[[col for col in csf_preprocessed.columns if 'VECPAC' in col]]
csf_lps = csf_preprocessed[[col for col in csf_preprocessed.columns if 'LPS' in col]]

csf_dss = rename_csf_columns(csf_dss, sara_ids)
csf_vecpac = rename_csf_columns(csf_vecpac, sara_ids)
csf_lps = rename_csf_columns(csf_lps, sara_ids)

common_dss = sorted(list(set(str_dss.columns) & set(csf_dss.columns)))
common_vecpac = sorted(list(set(str_vecpac.columns) & set(csf_vecpac.columns)))
common_lps = sorted(list(set(str_lps.columns) & set(csf_lps.columns)))

print(common_dss, common_vecpac, common_lps)

# simplify to the mice they have in common
str_aligned_dss = str_dss[["ID"] + common_dss]
str_aligned_vecpac = str_vecpac[["ID"] + common_vecpac]
str_aligned_lps = str_lps[["ID"] + common_lps]
csf_aligned_dss = csf_dss[common_dss]
csf_aligned_vecpac = csf_vecpac[common_vecpac]
csf_aligned_lps = csf_lps[common_lps]

csf_aligned_dss.insert(0, "ID", csf_aligned_dss.index)
csf_aligned_vecpac.insert(0, "ID", csf_aligned_vecpac.index)
csf_aligned_lps.insert(0, "ID", csf_aligned_lps.index)

# merge each model by mouse
merged_dss = pd.concat([str_aligned_dss, csf_aligned_dss], axis=0)
merged_vecpac = pd.concat([str_aligned_vecpac, csf_aligned_vecpac], axis=0)
merged_lps = pd.concat([str_aligned_lps, csf_aligned_lps], axis=0)

print("MERGED DSS: ", merged_dss.columns)
merged_dss.to_csv("merged_dss.csv", index=False)
print("MERGED VECPAC: ", merged_vecpac.columns)
merged_vecpac.to_csv("merged_vecpac.csv", index=False)
print("MERGED LPS: ", merged_lps.columns)
merged_lps.to_csv("merged_lps.csv", index=False)
