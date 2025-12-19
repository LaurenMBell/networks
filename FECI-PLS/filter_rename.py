import pandas as pd
import numpy as np
from scipy import stats

"""
What you're going to do:
1) Filter to control samples.
2) Combine PLS and CPX data tables (join by mouse).
3) Compute pooled correlations (for p-value, correlation coefficient)
4) Compute grouped correlation (for correlation coefficient) (within VECPAC, LPS, & DSS)
5) Filter to edges where sign(pooled) = sign(VECPAC) = sign(LPS) = sign(DSS)
6) FDR correct correlation p-values
"""

# tokenize, select the right ID (the short one) 
def rename_columns(df):
    new_cols = []
    prefixes = ("32", "33", "36", "9", "7", "10", "T0", "89") 

    for col in df.columns:
        tokens = col.split("_")
        match = next((t for t in tokens if t.startswith(prefixes)), None)

        if match:
            new_cols.append(match)
        else:
            new_cols.append(col)

    df.columns = new_cols
    return df


pls_dss = pd.read_csv("pls_data/pls_DSS_post_norm.csv")
pls_vecpac = pd.read_csv("pls_data/pls_VECPAC_post_norm.csv")
pls_lps = pd.read_csv("pls_data/pls_LPS_post_norm.csv")

feci_dss = pd.read_csv("feci_data/feci_DSS_post_norm.csv")
feci_vecpac = pd.read_csv("feci_data/feci_VECPAC_post_norm.csv")
feci_lps = pd.read_csv("feci_data/feci_LPS_post_norm.csv")

#============= TOKENIZE NAMES ===========================

#change the column names of all the mice
feci_dss = rename_columns(feci_dss)
feci_vecpac = rename_columns(feci_vecpac)
feci_lps = rename_columns(feci_lps)
pls_dss = rename_columns(pls_dss)
pls_vecpac = rename_columns(pls_vecpac)
pls_lps = rename_columns(pls_lps)

print("PLS - DSS:", pls_dss.columns)
print("FECI - DSS:",feci_dss.columns)
print("PLS - VECPAC:",pls_vecpac.columns)
print("FECI - VECPAC:",feci_vecpac.columns)
print("PLS - LPS:",pls_lps.columns)
print("FECI - LPS:", feci_lps.columns)

#get rid of the control row
pls_dss = pls_dss.drop(0)
pls_vecpac = pls_vecpac.drop(0)
pls_lps = pls_lps.drop(0)

feci_dss = feci_dss.drop(0)
feci_vecpac = feci_vecpac.drop(0)
feci_lps = feci_lps.drop(0)

# annotate feci metabolites
feci_dss["ID"]    = feci_dss["ID"].astype(str).str.strip() + "-F"
feci_vecpac["ID"] = feci_vecpac["ID"].astype(str).str.strip() + "-F"
feci_lps["ID"]    = feci_lps["ID"].astype(str).str.strip() + "-F"

# annotate pls metabolites
pls_dss["ID"]    = pls_dss["ID"].astype(str).str.strip() + "-P"
pls_vecpac["ID"] = pls_vecpac["ID"].astype(str).str.strip() + "-P"
pls_lps["ID"]    = pls_lps["ID"].astype(str).str.strip() + "-P"


#=============== JOIN BY MOUSE ==========================
common_dss = list(set(feci_dss.columns) & set(pls_dss.columns))
common_vecpac = list(set(feci_vecpac.columns) & set(pls_vecpac.columns))
common_lps = list(set(feci_lps.columns) & set(pls_lps.columns))

#simplify to the mice they have in common
feci_aligned_dss = feci_dss[common_dss]
feci_aligned_vecpac = feci_vecpac[common_vecpac]
feci_aligned_lps = feci_lps[common_lps]
pls_aligned_dss = pls_dss[common_dss]
pls_aligned_vecpac = pls_vecpac[common_vecpac]
pls_aligned_lps = pls_lps[common_lps]

#merge each model by mouse
merged_dss = pd.concat([feci_aligned_dss, pls_aligned_dss], axis=0)
merged_vecpac = pd.concat([feci_aligned_vecpac, pls_aligned_vecpac], axis=0)
merged_lps = pd.concat([feci_aligned_lps, pls_aligned_lps], axis=0)

#clean ID column 
merged_dss["ID"] = merged_dss["ID"].str.strip()
merged_vecpac["ID"] = merged_vecpac["ID"].str.strip()
merged_lps["ID"] = merged_lps["ID"].str.strip()

print("MERGED DSS: ", merged_dss.columns)
merged_dss.to_csv("merged_dss.csv", index=False)
print("MERGED VECPAC: ", merged_vecpac.columns)
merged_vecpac.to_csv("merged_vecpac.csv", index=False)
print("MERGED LPS: ", merged_lps.columns)
merged_lps.to_csv("merged_lps.csv", index=False)



