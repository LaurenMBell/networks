import pandas as pd
import pubchempy as pcp

mets = pd.read_excel("metabolites_for_sara.xlsx")

cids = []
for name in mets["Metabolite Name"]:
    results = pcp.get_compounds(name, "name")
    if results:
        print(f"{name}: {results[0].cid}")
        cids.append({"Metabolite Name": name, "PubChem": results[0].cid, "SMILES":results[0].smiles})
    else:
        #try: 
            #suffixes = name.strip("")
            #name = 
        print(f"SKIPPED: {name}")
        cids.append({"Metabolite Name": name, "PubChem": None, "SMILES":None})

results_df = pd.DataFrame(cids)
results_df.to_csv("metabolites_for_sara.csv", index=False)

