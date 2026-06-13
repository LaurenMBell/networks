# 5-4-26
# script for comparing Lauren de-derivatized metabolite names against Marcello names
# most of this code comes from pubchem_rest.py

import pandas as pd
import pubchempy as pcp
import time

data = pd.read_excel("metabolites_analysis.xlsx", sheet_name="Marcello Comparison")
M = data["Marcello Name"]
M.attrs["Name"] = "M"

L = data["Lauren Name"]
L.attrs["Name"] = "L"

to_check = [M, L]

for net in to_check:
    mets = []  
    for name in net:
        try:
            results = pcp.get_compounds(name, "name")
        except:
            time.sleep(3)
            results = pcp.get_compounds(name, "name")

        if results:
            pubchem_id = str(results[0].cid)
            synonyms = results[0].synonyms
            print(f"{name}: {pubchem_id}")
            mets.append({'metabolite': name, 'pubchem_id': pubchem_id})  
                    
        else:
            print(f"SKIPPED: {name}")
            mets.append({'metabolite': name, 'pubchem_id': None}) 

    results_df = pd.DataFrame(mets)
    results_df.to_csv(f"m_comp_{net.attrs["Name"]}.tsv", sep = "\t", index=False)