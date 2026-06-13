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
    pd.read_csv("VECPAC_correlations.csv"), "VECPAC")
lps = rename_correlation_columns(
    pd.read_csv("LPS_correlations.csv"), "LPS")
dss = rename_correlation_columns(
    pd.read_csv("DSS_correlations.csv"), "DSS")
pooled = rename_correlation_columns(
    pd.read_csv("Pooled_correlations.csv"), "Pooled")

dss_subset = dss[["DSS_r", "DSS_p", "DSS_n"]]
lps_subset = lps[["LPS_r", "LPS_p", "LPS_n"]]
pooled_subset = pooled[["Pooled_r", "Pooled_p", "Pooled_n"]]

df = pd.concat([vecpac, dss_subset, lps_subset, pooled_subset], axis=1)

model_r_cols = ["VECPAC_r", "DSS_r", "LPS_r"]
all_r_cols = model_r_cols + ["Pooled_r"]

df["r Count"] = df[model_r_cols].notna().sum(axis=1)
df.to_csv("all_correlations_pre_FDR.csv", index=False)
print(f"saved all_correlations_pre_FDR.csv")


is_consistent = df[df["Pooled_p"].notna()].copy()

is_consistent["Sign"] = np.where(is_consistent["Pooled_r"] > 0, "+", "-")

valid_pvals_mask = pd.to_numeric(is_consistent["Pooled_p"], errors='coerce').notna()
valid_pvals = is_consistent.loc[valid_pvals_mask, "Pooled_p"]

is_consistent["Pooled FDR"] = np.nan

if len(valid_pvals) > 0:
    rejected, corrected_pvals = fdrcorrection(
        valid_pvals.values,
        alpha=0.05,
        method="indep",
        is_sorted=False
    )
    is_consistent.loc[valid_pvals_mask, "Pooled FDR"] = corrected_pvals

is_consistent.to_csv("CSF-RBN_FDR.csv", index=False)
print(f"saved CSF-RBN_FDR.csv")

thresholded = is_consistent[is_consistent["Pooled FDR"] <= 0.05]
thresholded.to_csv("CSF-RBN_FDR_threshold.csv", index=False)