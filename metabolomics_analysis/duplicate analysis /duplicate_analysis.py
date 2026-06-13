# 5-4-26
# script for comparing de-derivatized metabolites against existing metabolite list 
# result was processed in excel for repeating values 
# most of this code comes from pubchem_rest.py

import pandas as pd
import pubchempy as pcp
import time

PLS = pd.read_excel("metabolites_analysis.xlsx", sheet_name="Duplicate Check PLS")
PLS.attrs["Name"] = "PLS"
FEC = pd.read_excel("metabolites_analysis.xlsx", sheet_name="Duplicate Check FEC")
FEC.attrs["Name"] = "FEC"
CSF = pd.read_excel("metabolites_analysis.xlsx", sheet_name="Duplicate Check CSF")
CSF.attrs["Name"] = "CSF"

#networks = [PLS, FEC, CSF]
networks = [CSF]

for net in networks:
    mets = []  
    for name in net.iloc[:,0]:
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
    results_df.to_csv(f"duplicate_analysis_{net.attrs["Name"]}.tsv", sep = "\t", index=False)