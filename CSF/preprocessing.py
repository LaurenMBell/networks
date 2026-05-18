"""
STEP 1: 
     - Raw data in Box, PROJECTS.Brain-gut_SaraCarloni.Metabolomics
     - Drop samples according to QC
     - Filter for CTRLs
     - Log transformation
     - Quantile normalization
     - Drop 0s for NaNs
"""

import pandas as pd
import numpy as np
import qnorm


def log_1(x):
    return np.log2(x + 1)


#TO KEEP
keep_samples = [
    "DSS_76_090224",
    "DSS_B896_090224",
    "DSS_T0.7_090224",
    "DSS_T0.8_090224",
    "DSS_73_090224",
    "935_CTR_LPS_CSF_070224",
    "936_CTR_LPS_CSF_070224_2",
    "937_CTR_LPS_CSF_070224",
    "941_CTR_LPS_CSF_070224",
    "942_CTR_LPS_CSF_070224",
    "943_CTR_LPS_CSF_070224",
    "3297_VECPAC_190224",
    "3679_VECPAC_190224",
    "3681_VECPAC_190224"]

df = pd.read_csv("data/CSF_raw.csv", header=None)

samples = df.iloc[0, 2:].str.strip().str.replace(r'\s+', '_', regex=True).values
metabolites = df.iloc[2:, 0].values
data = df.iloc[2:, 2:].values

data = pd.DataFrame(data, index=metabolites, columns=samples)
data = data.apply(pd.to_numeric, errors="coerce")

data = data[[col for col in data.columns if col in keep_samples]]

data.to_csv("pre_norm_CSF.csv", index=True)

log_transformed = data.apply(log_1)
print(log_transformed[0:5])

log_transformed.to_csv("log_trans_CSF.csv")


lg_DSS = log_transformed[["DSS_76_090224","DSS_B896_090224",
    "DSS_T0.7_090224","DSS_T0.8_090224", "DSS_73_090224"]]
q_DSS = qnorm.quantile_normalize(lg_DSS, axis=1)
DSS_nans = q_DSS.copy()
DSS_nans[data[["DSS_76_090224","DSS_B896_090224",
    "DSS_T0.7_090224","DSS_T0.8_090224", "DSS_73_090224"]] == 0] = np.nan
DSS_nans.to_csv("CSF_DSS_preprocessed.csv", index=True)


lg_LPS = log_transformed[["935_CTR_LPS_CSF_070224",
    "936_CTR_LPS_CSF_070224_2", "937_CTR_LPS_CSF_070224",
    "941_CTR_LPS_CSF_070224", "942_CTR_LPS_CSF_070224",
    "943_CTR_LPS_CSF_070224"]]
q_LPS = qnorm.quantile_normalize(lg_LPS, axis=1)
LPS_nans = q_LPS.copy()
LPS_nans[data[["935_CTR_LPS_CSF_070224",
    "936_CTR_LPS_CSF_070224_2", "937_CTR_LPS_CSF_070224",
    "941_CTR_LPS_CSF_070224", "942_CTR_LPS_CSF_070224",
    "943_CTR_LPS_CSF_070224"]] == 0] = np.nan
LPS_nans.to_csv("CSF_LPS_preprocessed.csv", index=True)


lg_VECPAC = log_transformed[["3297_VECPAC_190224",
    "3679_VECPAC_190224","3681_VECPAC_190224"]]
q_VECPAC = qnorm.quantile_normalize(lg_VECPAC, axis=1)
VECPAC_nans = q_VECPAC.copy()
VECPAC_nans[data[["3297_VECPAC_190224",
    "3679_VECPAC_190224","3681_VECPAC_190224"]] == 0] = np.nan
VECPAC_nans.to_csv("CSF_VECPAC_preprocessed.csv", index=True)

