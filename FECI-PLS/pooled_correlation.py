import pandas as pd
from scipy import stats
import itertools
import numpy as np
import chime

def compute_correlations(data, model_name):
    nodes = data['ID'].values

    nodes_clean = (
    pd.Series(nodes)
      .astype(str)
      .str.strip()         
      .str.replace("–", "-", regex=False)
      .values
    )
    
    #identify genes and metabolites
    p_data = [
        (i, node_val)
        for i, node_val in enumerate(nodes_clean)
        if node_val.endswith("-P")
    ]

    f_data = [
        (i, node_val)
        for i, node_val in enumerate(nodes_clean)
        if node_val.endswith("-F")
    ]

    
    #get sample columns 
    mice = [col for col in data.columns if col != "ID"]
    
    #create all metabolite-gene pairs
    pairs = list(itertools.product(p_data, f_data))
        
    #pairs = ['(2,2-Dimethyl-5-[2-(2-ethoxymethoxy)propyl][1,3]dioxolan-4-yl)methanol-F'), (1,'(7a-Isopropenyl-4,5-dimethyloctahydroinden-4-yl)methanol-P')]

    results = []
    
    #for (p_idx, p_id), (f_idx, f_id) in pairs:
    for (p_idx, p_id), (f_idx, f_id) in pairs:
        p_values = data.iloc[p_idx][mice].astype(float).to_numpy()
        f_values = data.iloc[f_idx][mice].astype(float).to_numpy()
        
        mask = ~np.isnan(p_values) & ~np.isnan(f_values)
        p_masked = p_values[mask]
        f_masked = f_values[mask]
        n = mask.sum()


        if n < 3: #if sample size is less than 3
            corr, pval = np.nan, np.nan 
        elif np.all(p_masked == p_masked[0]) or np.all(f_masked == f_masked[0]):
            corr, pval = np.nan, np.nan
        else: 
            #print(f"{c_id} - {f_id}\n")

            corr, pval = stats.spearmanr(p_masked, f_masked, nan_policy='omit')
            
            if abs(corr) == 1:
                pval = 0
        
        results.append({
            'pls_metabolite': p_id, 
            'feci_metabolite': f_id, 
            f'{model_name}_corr_coef': corr, 
            f'{model_name}_p-value': pval, 
            f'{model_name}_n_samples': n
        })
    
    return pd.DataFrame(results)

#compute correlations for each model separately
dss = pd.read_csv("merged_dss.csv")
dss_correlations = compute_correlations(dss, "DSS")
dss_correlations.to_csv("DSS_correlations.csv", index=False)
print("dss done")
chime.success()

lps = pd.read_csv("merged_lps.csv")
lps_correlations = compute_correlations(lps, "LPS")
lps_correlations.to_csv("LPS_correlations.csv", index=False)
print("lps done")
chime.success()


vecpac = pd.read_csv("merged_vecpac.csv")
vecpac_correlations = compute_correlations(vecpac, "VECPAC")
vecpac_correlations.to_csv("VECPAC_correlations.csv", index=False)
print("vecpac done")
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
print("pooled done")
chime.success()