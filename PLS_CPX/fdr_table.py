import scipy.stats
import numpy as np
import pandas as pd
from statsmodels.stats.multitest import fdrcorrection
import chime

# Load CSVs
vecpac = pd.read_csv("VECPAC_correlations.csv")
lps    = pd.read_csv("LPS_correlations.csv")
dss    = pd.read_csv("DSS_correlations.csv")
pooled = pd.read_csv("pooled_correlations.csv")

#print(lps)
# Rename correlation & p-value columns
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

# Select only DSS r/p-values and LPS r/p-values
dss_subset = dss[['DSS r', 'DSS p-values', 'DSS n']]
lps_subset = lps[['LPS r', 'LPS p-values', 'LPS n']]
pooled_subset = pooled[['Pooled r', 'Pooled p-values', 'Pooled n']]

# Concatenate columns and drop Nan values 
df = pd.concat([vecpac, dss_subset, lps_subset, pooled_subset], axis=1).dropna()

# Add a new column based on sign consistency
df['Consistent'] = (
    ((df['VECPAC r'] > 0) & (df['LPS r'] > 0) & (df['DSS r'] > 0) & (df['Pooled r'] > 0)) |
    ((df['VECPAC r'] < 0) & (df['LPS r'] < 0) & (df['DSS r'] < 0) & (df['Pooled r'] < 0))
)

is_consistent = df.loc[df["Consistent"] == True].copy()

is_consistent["Sign"] = np.where(is_consistent['Pooled r'] > 0, "+", "-")

rejected, corrected_pvals = fdrcorrection(is_consistent.loc[:, "Pooled p-values"], alpha=0.05, method='indep',is_sorted=False)

is_consistent["Pooled FDR"] = corrected_pvals

is_consistent.to_csv("PLS-CPX_FDR.csv", index=False)
chime.success()