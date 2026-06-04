import numpy as np
import pandas as pd
from statsmodels.stats.multitest import fdrcorrection


def rename_correlation_columns(df, model_name):
    df = df.copy()
    df.columns.values[2] = f"{model_name}_r"
    df.columns.values[3] = f"{model_name}_p"
    df.columns.values[4] = f"{model_name}_n"
    return df


vecpac = rename_correlation_columns(
    pd.read_csv("correlations/VECPAC_correlations.csv"), "VECPAC")
lps = rename_correlation_columns(
    pd.read_csv("correlations/LPS_correlations.csv"), "LPS")
dss = rename_correlation_columns(
    pd.read_csv("correlations/DSS_correlations.csv"), "DSS")
pooled = rename_correlation_columns(
    pd.read_csv("correlations/Pooled_correlations.csv"), "Pooled")

dss_subset = dss[["DSS_r", "DSS_p", "DSS_n"]]
lps_subset = lps[["LPS_r", "LPS_p", "LPS_n"]]
pooled_subset = pooled[["Pooled_r", "Pooled_p", "Pooled_n"]]

df = pd.concat([vecpac, dss_subset, lps_subset, pooled_subset], axis=1)

model_r_cols = ["VECPAC_r", "DSS_r", "LPS_r"]
all_r_cols = model_r_cols + ["Pooled_r"]

df["r Count"] = df[model_r_cols].notna().sum(axis=1)
df_valid = df[df["r Count"] >= 2].copy()

df_valid["Consistent"] = df_valid[all_r_cols].apply(
    lambda row: (np.nanmin(row) > 0) or (np.nanmax(row) < 0), axis=1)

df_valid.to_csv("all_correlations_pre_FDR.csv", index=False)
print("saved all_correlations_pre_FDR.csv")

is_consistent = df_valid[
    df_valid["Consistent"] & df_valid["Pooled_p"].notna()].copy()

is_consistent["Sign"] = np.where(is_consistent["Pooled_r"] > 0, "+", "-")

i, corrected_pvals = fdrcorrection(is_consistent["Pooled_p"], alpha=0.05, method="indep", is_sorted=False)
is_consistent["Pooled FDR"] = corrected_pvals

is_consistent.to_csv("CSF_FDR.csv", index=False)
print("saved CSF_FDR.csv")

thresholded = is_consistent[is_consistent["Pooled FDR"] <= 0.05]
thresholded.to_csv("CSF_FDR_threshold.csv", index=False)