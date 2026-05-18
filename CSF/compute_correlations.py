"""
STEP 2: Compute all correlations between CSF-CSF
"""

import pandas as pd
import numpy as np
import itertools
from scipy import stats


def correlations(data, label):
    metabolites = data.index.tolist()
    pairs = list(itertools.combinations(metabolites, 2))
    results = []
    
    for met1, met2 in pairs:
        row1 = data.loc[met1].values.astype(float)
        row2 = data.loc[met2].values.astype(float)
        
        #NaNs
        df = pd.DataFrame({'x': row1, 'y': row2}).dropna()
        
        n = len(df)
        if n < 3:
            results.append((met1, met2, np.nan, np.nan, 0))
        else:
            corr, pval = stats.spearmanr(df['x'], df['y'])
            results.append((met1, met2, corr, pval, n))
    
    return results

"""
data = pd.read_csv("CSF_preprocessed.csv", index_col=0)

#split by model name
dss_cols = [col for col in data.columns if "DSS" in col]
vecpac_cols = [col for col in data.columns if "VECPAC" in col]
lps_cols = [col for col in data.columns if "LPS" in col]

dss_data = data[dss_cols]
vecpac_data = data[vecpac_cols]
lps_data = data[lps_cols]

dss_data.to_csv("data/CSF_DSS.csv")
vecpac_data.to_csv("data/CSF_VECPAC.csv")
lps_data.to_csv("data/CSF_LPS.csv")
"""
dss_data = pd.read_csv("CSF_DSS_preprocessed.csv", index_col=0)
vecpac_data = pd.read_csv("CSF_VECPAC_preprocessed.csv", index_col=0)
lps_data = pd.read_csv("CSF_LPS_preprocessed.csv", index_col=0)
data = pd.concat([dss_data, vecpac_data, lps_data], axis=1)


models = [("DSS", dss_data),("VECPAC", vecpac_data),
    ("LPS", lps_data),("Pooled", data)]

for model_name, model_data in models:
    results = correlations(model_data, model_name)

    df_results = pd.DataFrame(results, columns=["n1", "n2", "r", "p-value", "n"])
    df_results.to_csv(f"correlations/{model_name}_correlations.csv", index=False)
    
    print(f"done with {model_name}")