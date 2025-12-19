import scipy.stats
import numpy as np
import pandas as pd
from statsmodels.stats.multitest import fdrcorrection
import chime

#load CSVs
vecpac = pd.read_csv("VECPAC_correlations.csv")
lps    = pd.read_csv("LPS_correlations.csv")
dss    = pd.read_csv("DSS_correlations.csv")
pooled = pd.read_csv("pooled_correlations.csv")

# get rid of the tissue flags for the edge table
for df in [vecpac, lps, dss, pooled]: 
    df.iloc[:, 0] = df.iloc[:, 0].astype(str).str.replace(r"-P$", "", regex=True)
    df.iloc[:, 1] = df.iloc[:, 1].astype(str).str.replace(r"-F$", "", regex=True)

#print(lps)
#rename correlation & p-value columns
vecpac.columns.values[2] = 'VECPAC r'
vecpac.columns.values[3] = 'VECPAC p-values'
vecpac.columns.values[4] = 'VECPAC n'
lps.columns.values[2]    = 'LPS r'
lps.columns.values[3]    = 'LPS p-values'
lps.columns.values[4] = 'LPS n'
dss.columns.values[2]    = 'DSS r'
dss.columns.values[3]    = 'DSS p-values'
dss.columns.values[4] = 'DSS n'
pooled.columns.values[2] = 'Pooled r'
pooled.columns.values[3] = 'Pooled p-values'
pooled.columns.values[4] = 'Pooled n'

#select only DSS r/p-values and LPS r/p-values
dss_subset = dss[['DSS r', 'DSS p-values', 'DSS n']]
lps_subset = lps[['LPS r', 'LPS p-values', 'LPS n']]
pooled_subset = pooled[['Pooled r', 'Pooled p-values', 'Pooled n']]

#concatenate columns and drop Nan values 
df = pd.concat([vecpac, dss_subset, lps_subset, pooled_subset], axis=1)

r_cols = ["VECPAC r", "DSS r", "LPS r"]
df["r Count"] = df[r_cols].notna().sum(axis=1)

# require >= 2 models to be non-nan
df_valid = df[df["r Count"] >= 2].copy()

r_cols = ["VECPAC r", "DSS r", "LPS r", "Pooled r"]

df_valid["Consistent"] = df_valid[r_cols].apply(
    lambda row: (np.nanmin(row) > 0) or (np.nanmax(row) < 0), axis=1
)

df_valid.to_csv("all_correlations_pre_FDR.csv", index=False)
print("saved all_correlations_pre_FDR.csv")

is_consistent = df_valid[df_valid["Consistent"]].copy()

is_consistent["Sign"] = np.where(is_consistent['Pooled r'] > 0, "+", "-")

rejected, corrected_pvals = fdrcorrection(is_consistent.loc[:, "Pooled p-values"], alpha=0.05, method='indep',is_sorted=False)

is_consistent["Pooled FDR"] = corrected_pvals

is_consistent.to_csv("FECI-PLS_FDR.csv", index=False)
chime.success()