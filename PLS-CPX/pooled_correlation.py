import pandas as pd
from scipy import stats
import itertools
import numpy as np
import chime

def compute_correlations(data, model_name):
    nodes = data['ID'].values
    
    #identify genes and metabolites
    gene_data = [(i, node_val) for i, node_val in enumerate(nodes) 
                 if (str(node_val).startswith("ENSMUSG") and str(node_val).endswith("-P"))]
    metabolite_data = [(i, node_val) for i, node_val in enumerate(nodes) 
                       if not str(node_val).startswith("ENSMUSG")]
    
    #get sample columns 
    mice = [col for col in data.columns if col != "ID"]
    
    #create all metabolite-gene pairs
    pairs = list(itertools.product(metabolite_data, gene_data))
        
    results = []
    
    for (p_idx, p_id), (c_idx, c_id) in pairs:
        p_values = data.iloc[p_idx][mice].astype(float).to_numpy()
        c_values = data.iloc[c_idx][mice].astype(float).to_numpy()
        
        mask = ~np.isnan(p_values) & ~np.isnan(c_values)
        p_masked = p_values[mask]
        c_masked = c_values[mask]
        n = mask.sum()


        if n < 3: #if sample size is less than 3
            corr, pval = np.nan, np.nan 
        elif np.all(p_masked == p_masked[0]) or np.all(c_masked == c_masked[0]):
            corr, pval = np.nan, np.nan
        else: 
            #print(f"{c_id} - {f_id}\n")

            corr, pval = stats.spearmanr(p_masked, c_masked, nan_policy='omit')
            
            if abs(corr) == 1:
                pval = 0
        
        results.append({
            'cpx_gene': c_id, 
            'pls_metabolite': p_id,
            f'{model_name}_corr_coef': corr, 
            f'{model_name}_p-value': pval, 
            f'{model_name}_n_samples': n
        })
    
    return pd.DataFrame(results)

#read data
dss = pd.read_csv("merged_dss.csv")
dss_correlations = compute_correlations(dss, "DSS")
dss_correlations.to_csv("DSS_correlations.csv", index=False)
print("dss done")
chime.success()

vecpac = pd.read_csv("merged_vecpac.csv")
vecpac_correlations = compute_correlations(vecpac, "VECPAC")
vecpac_correlations.to_csv("VECPAC_correlations.csv", index=False)
print("vecpac done")
chime.success()

lps = pd.read_csv("merged_lps.csv")
lps_correlations = compute_correlations(lps, "LPS")
lps_correlations.to_csv("LPS_correlations.csv", index=False)
print("lps done")
chime.success()

# ============= POOLED CORRELATIONS =============#


#rename columns to distinguish samples by model
dss_renamed = dss.rename(columns={col: f"{col}_DSS" if col != 'ID' else col for col in dss.columns})
vecpac_renamed = vecpac.rename(columns={col: f"{col}_VECPAC" if col != 'ID' else col for col in vecpac.columns})
lps_renamed = lps.rename(columns={col: f"{col}_LPS" if col != 'ID' else col for col in lps.columns})

#drop extra ID columns
vecpac_renamed.drop("ID", axis=1, inplace=True)
lps_renamed.drop("ID", axis=1, inplace=True)

#concatenate all models
all_data = pd.concat([dss_renamed, vecpac_renamed, lps_renamed], axis=1)

all_data.to_csv("merged_all_data.csv", index=False)

#compute pooled correlationss
pooled_correlations = compute_correlations(all_data, "pooled")
pooled_correlations.to_csv("pooled_correlations.csv", index=False)
chime.success()