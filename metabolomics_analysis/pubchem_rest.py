# late April/early May 2026
# script to pull PubChem and HMDB IDs for metabolites 
#takes blank workbork to fill out, looks up metablite in pubchem, 

import pandas as pd
import pubchempy as pcp

mets = pd.read_excel("metabolites_for_sara.xlsx")

mapping = pd.read_csv("pubchemID_mapping.csv", dtype=str)
mapping = mapping.set_index("Query")

cids = []
for name in mets["Metabolite Name"]:
    results = pcp.get_compounds(name, "name")
    if results:
        pubchem_id = str(results[0].cid)
        hmdb_code = None
        if pubchem_id in mapping.index:
            hmdb_code = mapping.at[pubchem_id, "HMDB"]
        print(f"{name}: {pubchem_id}")
        cids.append({
            "Metabolite Name": name, "PubChem": pubchem_id,
            "HMDB": hmdb_code, "SMILES": results[0].smiles})
    else:
        print(f"SKIPPED: {name}")
        cids.append({"Metabolite Name": name, "PubChem": None,
            "HMDB": None,"SMILES": None})

results_df = pd.DataFrame(cids)
results_df.to_csv("metabolites_for_sara.csv", index=False)

