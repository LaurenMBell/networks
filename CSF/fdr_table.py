import numpy as np
import scipy
import pandas as pd
from statsmodels.stats.multitest import fdrcorrection


def metacor(r, n, do_fisher=True, return_transformed=True, cor_type="pearson", meta_type="fixed"):
	if cor_type not in ["pearson", "spearman"]:
		raise Exception(f"Invalid correlation type \"{cor_type}\"")
	if meta_type not in ["fixed", "random"]:
		raise Exception(f"Invalid meta analysis type \"{meta_type}\"")
	if not do_fisher:
		raise Exception("Non-Fisher's-z-transformed correlations not yet supported")
	if meta_type == "random":
		raise Exception("Random effects meta analysis not yet supported")
	if not return_transformed:
		raise Exception("Returning an untransformed common correlation coefficient estimate not yet supported")

	do_override = False
	if do_fisher:
		# Explicitly handle cases of r = 1 or r = -1
		r_p_one = (r == 1.0).any(axis=1)
		r_n_one = (r == -1.0).any(axis=1)
		override = (r_p_one | r_n_one)
		do_override = np.any(override)
		if do_override:
			# r_common = inf or -inf depending on sign of r
			r_common_override_val = np.where(r_p_one, np.inf, np.where(r_n_one, -np.inf, np.nan))
			# p = 0
			p_override_val = np.where(override, 0.0, np.nan)
			# Convert `override` to a column of a matrix, then broadcast to match `r`
			override_matrix = np.broadcast_to(np.reshape(override, (-1, 1)), r.shape)
			r = np.where(override_matrix, np.nan, r)
			n = np.where(override_matrix, np.nan, n)
		r = np.arctanh(r)

	if cor_type == "pearson":
		if do_fisher:
			# Bonett (2000)
			variance = 1 / (n - 3)
		else:
			raise Exception("NYI")
	elif cor_type == "spearman":
		if do_fisher:
			# Bonett (2000)
			variance = (1 + (np.square(r) / 2)) / (n - 3)
		else:
			raise Exception("NYI")
	# Cooper (2009)
	weight = 1 / variance
	variance_common = 1 / np.nansum(weight, axis=1)

	r_common = np.nansum(np.multiply(weight, r), axis=1) * variance_common
	stddev_common = np.sqrt(variance_common)

	z = r_common / stddev_common
	p = 2 * scipy.stats.norm.sf(np.abs(z))

	if do_override:
		r_common = np.where(override, r_common_override_val, r_common)
		p = np.where(override, p_override_val, p)

	return r_common, p

# Load CSVs
vecpac = pd.read_csv("correlations/VECPAC_correlations.csv")
lps    = pd.read_csv("correlations/LPS_correlations.csv")
dss    = pd.read_csv("correlations/DSS_correlations.csv")

# Rename correlation & p-value columns
vecpac.columns.values[2] = 'VECPAC_r'
vecpac.columns.values[3] = 'VECPAC_p'
vecpac.columns.values[4] = "VECPAC_n"
lps.columns.values[2]    = 'LPS_r'
lps.columns.values[3]    = 'LPS_p'
lps.columns.values[4] = "LPS_n"
dss.columns.values[2]    = 'DSS_r'
dss.columns.values[3]    = 'DSS_p'
dss.columns.values[4] = "DSS_n"


# Select only DSS r/p-values and LPS r/p-values
dss_subset = dss[['DSS_r', 'DSS_p', "DSS_n"]]
lps_subset = lps[['LPS_r', 'LPS_p', 'LPS_n']]

#df = vecpac.merge(dss, on=["Metabolite 1", "Metabolite 2"], how="right")
#df = df.merge(lps, on=["Metabolite 1", "Metabolite 2"], how="right")


# Concatenate columns
df = pd.concat([vecpac, dss_subset, lps_subset], axis=1)
#df = pd.concat([vecpac, dss, lps], axis=1)

#df = (
    #vecpac
    #.merge(dss, left_index=False, right_index=False, left_on=["Metabolite 1", "Metabolite 2"], right_on=["Metabolite 1", "Metabolite 2"], how="outer")
    #.merge(lps, on=["Metabolite 1", "Metabolite 2"], how="outer")
#)

# Add a new column based on sign consistency
df['Consistent'] = (
    ((df['VECPAC_r'] > 0) & (df['LPS_r'] > 0) & (df['DSS_r'] > 0)) |
    ((df['VECPAC_r'] < 0) & (df['LPS_r'] < 0) & (df['DSS_r'] < 0))
)

to_check = ["VECPAC_r", "LPS_r", "DSS_r"]

#i give up on coming up with a cleverer solution ngl

df['r_values_count'] = df[to_check].notna().sum(axis=1)

df["Meta-analysis_validity"] = ((df[to_check].notna().sum(axis=1) >= 2) &
    (np.sign(df[to_check]).sum(axis=1).abs() == df[to_check].notna().sum(axis=1))
)
    

# Filter rows where Consistent? == True
is_valid = df["Meta-analysis_validity"]
df_valid = df.loc[is_valid].copy()

# get the correlation coefficients and respective array sizes for metacor
r = df_valid[['VECPAC_r', 'DSS_r', 'LPS_r']].to_numpy()
n = df_valid[['VECPAC_n', 'DSS_n', 'LPS_n']].to_numpy()

# Run meta-analysis
r_common, p_values = metacor(r, n, do_fisher=True, cor_type="spearman", meta_type="fixed")

# Add raw p-values column only to consistent metabolites
df_valid['meta_p'] = p_values
# Add meta-analysis values into original table
df['meta_p'] = np.nan
df.loc[is_valid, 'meta_p'] = p_values

# Run FDR correction only on valid p-values
is_valid = ~df['meta_p'].isna()
rejected, pvals_corrected = fdrcorrection(df.loc[is_valid, 'meta_p'], alpha=0.05, method='indep',is_sorted=False)

# Add FDR results
df['FDR'] = np.nan
df.loc[is_valid, 'FDR'] = pvals_corrected
df = df.drop('r_values_count', axis =1)

# Save final file with meta-analysis and FDR correction
fdr_file = "CSF_FDR.csv"
df.to_csv(fdr_file, index=False)
print("saved FDR", fdr_file)



