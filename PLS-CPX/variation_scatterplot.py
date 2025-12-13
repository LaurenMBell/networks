import matplotlib.pyplot as plt
import pandas as pd 
import numpy as np

lps = pd.read_csv("metabolomics_data/pls_LPS_post_norm.csv")
vecpac = pd.read_csv("metabolomics_data/pls_VECPAC_post_norm.csv")
dss = pd.read_csv("metabolomics_data/pls_DSS_post_norm.csv")

models = [("LPS", lps), ("VECPAC", vecpac), ("DSS", dss)]

for name, model in models:
    model = model.drop(model.index[0]).reset_index(drop=True)
    
    for col in model.columns:
        if col not in ['ID', 'GROUP']:
            model[col] = pd.to_numeric(model[col], errors='coerce')

    samples = model.select_dtypes(include=[np.number]).columns.tolist()
    model_data = []

    for idx, row in model.iterrows():

        m_data = row[samples]
        
        if m_data.isna().all():
            metabolite_name = row.get('ID', f'Row {idx}')
            print(f"{metabolite_name} all NaN")
            continue
        

        m_mean = m_data.mean()
        deviations = m_data - m_mean
        m_median = np.abs(deviations.median())
        
        model_data.append(m_median)
    
    plt.figure()
    plt.hist(model_data, bins=50, color='blue', edgecolor='black')
    plt.xlabel("Metabolite Absolute Median Deviation")
    plt.ylabel("Counts")
    plt.title(f"{name} - Metabolite Deviation from Mean")
    plt.savefig(f"{name}_Median_Deviation.png")