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
    c_data = [
        (i, node_val)
        for i, node_val in enumerate(nodes_clean)
        if node_val.startswith("ENSMUSG") and node_val.endswith("-P")
    ]

    f_data = [
        (i, node_val)
        for i, node_val in enumerate(nodes_clean)
        if not node_val.startswith("ENSMUSG")
    ]

    
    #get sample columns 
    mice = [col for col in data.columns if col != "ID"]
    
    #create all metabolite-gene pairs
    pairs = list(itertools.product(c_data, f_data))
        
    results = []
    
    for (c_idx, c_id), (f_idx, f_id) in pairs:
        c_values = data.iloc[c_idx][mice].astype(float)
        f_values = data.iloc[f_idx][mice].astype(float)
        
        df = pd.DataFrame({'cpx_gene': c_values, 'feci_metabolite': f_values})
        n = (~(df["cpx_gene"].isna() | df["feci_metabolite"].isna())).sum()
        
        if n < 3:
            corr, pval = np.nan, np.nan  
        else:
            print(f"{c_id} - {f_id}\n")
            corr, pval = stats.spearmanr(df['cpx_gene'], df['feci_metabolite'], nan_policy='omit')
            
            if corr == 1 or corr == -1:
                pval = 0
        
        results.append({
            'cpx_gene': c_id, 
            'feci_metabolite': f_id, 
            f'{model_name}_corr_coef': corr, 
            f'{model_name}_p-value': pval, 
            f'{model_name}_n_samples': n
        })
    
    return pd.DataFrame(results)

#read data
dss = pd.read_csv("merged_dss.csv")
vecpac = pd.read_csv("merged_vecpac.csv")
lps = pd.read_csv("merged_lps.csv")

#compute correlations for each model separately
dss_correlations = compute_correlations(dss, "DSS")
dss_correlations.to_csv("DSS_correlations.csv", index=False)
chime.success()

vecpac_correlations = compute_correlations(vecpac, "VECPAC")
vecpac_correlations.to_csv("VECPAC_correlations.csv", index=False)
chime.success()

lps_correlations = compute_correlations(lps, "LPS")
lps_correlations.to_csv("LPS_correlations.csv", index=False)
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