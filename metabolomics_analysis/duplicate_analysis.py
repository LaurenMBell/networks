# 5-31-26
# script for comparing de-derivatized metabolites against existing metabolite list 
# result was processed in excel for repeating values 
# most of this code comes from pubchem_rest.py

import pandas as pd
import pubchempy as pcp

PLS = pd.read_excel("metabolites_analysis.xlsx", sheet_name="Duplicate Check PLS")
PLS.attrs["Name"] = "PLS"
FEC = pd.read_excel("metabolites_analysis.xlsx", sheet_name="Duplicate Check FEC")
FEC.attrs["Name"] = "FEC"
CSF = pd.read_excel("metabolites_analysis.xlsx", sheet_name="Duplicate Check CSF")
CSF.attrs["Name"] = "CSF"

networks = [PLS, FEC, CSF]

mets = {}
for net in networks:
    for name in net.iloc[:,0]:
        results = pcp.get_compounds(name, "name")
        if results:
            pubchem_id = str(results[0].cid)
            synonyms = results[0].synonyms
            print(f"{name}: {pubchem_id}")
            mets[name] = synonyms

            for met in mets:
                for syn in mets[synonyms]:
                    if name == syn:
                        print(f"DUPLICATE: {name} [{pubchem_id}]")
                    
        else:
            print(f"SKIPPED: {name}")
            mets[name] = None

    results_df = pd.DataFrame(mets)
    results_df.to_csv(f"duplicate_analysis_{net.attrs["Name"]}.csv", sep = "\t", index=False)