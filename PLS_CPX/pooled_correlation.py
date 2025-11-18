import pandas as pd
from scipy import stats
import itertools
import numpy as np
import chime

def compute_correlations(data, model_name):
    """
    Compute correlations for a single model (DSS, VECPAC, or LPS)
    """
    nodes = data['ID'].values
    
    # Identify genes and metabolites
    gene_data = [(i, node_val) for i, node_val in enumerate(nodes) 
                 if (str(node_val).startswith("ENSMUSG") and str(node_val).endswith("-P"))]
    metabolite_data = [(i, node_val) for i, node_val in enumerate(nodes) 
                       if not str(node_val).startswith("ENSMUSG")]
    
    # Get sample columns (exclude ID)
    mice = [col for col in data.columns if col != "ID"]
    
    # Create all metabolite-gene pairs
    pairs = list(itertools.product(metabolite_data, gene_data))
    
    print(f"{model_name}: Computing {len(pairs)} correlations with {len(mice)} samples")
    
    results = []
    
    for (met_idx, met_id), (gene_idx, gene_id) in pairs:
        met_values = data.iloc[met_idx][mice].astype(float)
        gene_values = data.iloc[gene_idx][mice].astype(float)
        
        df = pd.DataFrame({'metabolite': met_values, 'gene': gene_values}).dropna()
        n = len(df)
        
        if n < 3:
            corr, pval = np.nan, np.nan  
        else:
            corr, pval = stats.spearmanr(df['metabolite'], df['gene'])
        
        results.append({
            'metabolite': met_id, 
            'gene': gene_id, 
            f'{model_name}_corr_coef': corr, 
            f'{model_name}_p-value': pval, 
            f'{model_name}_n_samples': n
        })
    
    return pd.DataFrame(results)

# Read data
dss = pd.read_csv("merged_dss.csv")
vecpac = pd.read_csv("merged_vecpac.csv")
lps = pd.read_csv("merged_lps.csv")

# Compute correlations for each model separately
dss_correlations = compute_correlations(dss, "DSS")
vecpac_correlations = compute_correlations(vecpac, "VECPAC")
lps_correlations = compute_correlations(lps, "LPS")

# Save individual model results
dss_correlations.to_csv("DSS_correlations.csv", index=False)
vecpac_correlations.to_csv("VECPAC_correlations.csv", index=False)
lps_correlations.to_csv("LPS_correlations.csv", index=False)
chime.success()

# ============= POOLED CORRELATIONS =============

# Rename columns to distinguish samples by model
dss_renamed = dss.rename(columns={col: f"{col}_DSS" if col != 'ID' else col for col in dss.columns})
vecpac_renamed = vecpac.rename(columns={col: f"{col}_VECPAC" if col != 'ID' else col for col in vecpac.columns})
lps_renamed = lps.rename(columns={col: f"{col}_LPS" if col != 'ID' else col for col in lps.columns})

# Drop extra ID columns
vecpac_renamed.drop("ID", axis=1, inplace=True)
lps_renamed.drop("ID", axis=1, inplace=True)

# Concatenate all models
all_data = pd.concat([dss_renamed, vecpac_renamed, lps_renamed], axis=1)

# Compute pooled correlations
pooled_correlations = compute_correlations(all_data, "pooled")
pooled_correlations.to_csv("pooled_correlations.csv", index=False)
print(f"DSS: {len(dss_correlations)} correlations")
print(f"VECPAC: {len(vecpac_correlations)} correlations")
print(f"LPS: {len(lps_correlations)} correlations")
print(f"Pooled: {len(pooled_correlations)} correlations")
chime.success()