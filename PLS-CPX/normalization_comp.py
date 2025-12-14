import matplotlib.pyplot as plt
import pandas as pd 
import numpy as np

lps = pd.read_csv("metabolomics_data/pls_LPS_post_norm.csv")
vecpac = pd.read_csv("metabolomics_data/pls_VECPAC_post_norm.csv")
dss = pd.read_csv("metabolomics_data/pls_DSS_post_norm.csv")

#(name, model)
models = [("LPS", lps), ("VECPAC", vecpac), ("DSS", dss)]

for name, model in models:
    model = model.drop(model.index[0]).reset_index(drop=True)
    
    for col in model.columns:
        if col not in ['ID', 'GROUP']:
            model[col] = pd.to_numeric(model[col], errors='coerce')

    samples = model.select_dtypes(include=[np.number]).columns.tolist()
    #model_data = []
    df = pd.DataFrame()

    for idx, met in model.iterrows():

        m_data = met[samples]
        
        if m_data.isna().all():
            metabolite_name = met.get('ID', f'Row {idx}')
            print(f"{metabolite_name} all NaN")
            continue
        
            
        
        