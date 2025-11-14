import pandas as pd
import numpy as np
import re

"""
What you're going to do:
1) Filter to control samples.
2) Combine PLS and CPX data tables (join by mouse).
3) Compute pooled correlations (for p-value, correlation coefficient)
4) Compute grouped correlation (for correlation coefficient) (within VECPAC, LPS, & DSS)
5) Filter to edges where sign(pooled) = sign(VECPAC) = sign(LPS) = sign(DSS)
6) FDR correct correlation p-values
"""

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

pls_dss = pd.read_csv("metabolomics_data/pls_DSS_CTRL_post_norm.csv")
pls_vecpac = pd.read_csv("metabolomics_data/pls_VECPAC_CTRL_post_norm.csv")
pls_lps = pd.read_csv("metabolomics_data/pls_LPS_CTRL_post_norm.csv")

#============= CHANGE NAMES OF BOTH DATASETS ===========================

mouse_map = pd.read_csv("mouse_map.csv").to_dict()

#change the column names of the cpx mice
for df in [cpx_dss, cpx_vecpac, cpx_lps]:
     df.columns = [mouse_map.get(str(col), col) for col in df.columns]

#change the column names of the pls mive
pls_mice = set(pls_dss.columns) | set(pls_vecpac.columns) | set(pls_lps.columns) #set of all samples

for df in [pls_dss, pls_lps, pls_vecpac]:
    df.columns = [mouse_map.get(str(col), col) for col in df.columns]
   











