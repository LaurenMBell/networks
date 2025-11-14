import pandas as pd
import numpy as np
import qnorm

def log_1(x):
    return np.log2(x + 1)

def log_transform(infile, name):
    og_data = pd.read_csv(infile, header=[0,1])
    names = og_data.iloc[:,0]

    data_df = og_data.iloc[0:,1:]
    log_transformed = data_df.apply(log_1)

    named = log_transformed.copy()
    named.insert(0, "Names", names)
    named.to_csv(f"normalization/{name}_log_transformed.csv", index=False)

    #log_transformed.iloc[206,0:7] = np.nan
    #log_tran

    qnormed = qnorm.quantile_normalize(log_transformed, axis=1)
    qnormed.insert(0, "", names)
    qnormed.to_csv(f"normalization/{name}_qnormed.csv", index=False)

    qnormed_nans = qnormed.copy()
    qnormed_nans[log_transformed==0] = np.nan
    qnormed_nans.to_csv(f"normalization/{name}_post_norm.csv", index=False) #, na_rep="NaN")

log_transform("2. Filter_CTRL/feci_DSS_CTRL_filtered.csv", "feci_DSS")
log_transform("2. Filter_CTRL/feci_LPS_CTRL_filtered.csv", "feci_LPS")
log_transform("2. Filter_CTRL/feci_VECPAC_CTRL_filtered.csv", "feci_VECPAC")





