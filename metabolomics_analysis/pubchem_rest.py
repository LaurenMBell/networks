import pandas as pd
import pubchempy as pcp

mets = pd.read_excel("gut_brain_network_2026_02_17.xlsx", sheet_name="Metabolite External IDs")

records = []
for name in mets["Name"]:
    results = pcp.get_compounds(name, "name")
    if results:
        print(f"{name}: {results[0].cid}")
        records.append({"Name": name, "CID": results[0].cid})
    else:
        try: 
            suffixes = name.strip("")
            name = 
        print(f"SKIPPED: {name}")
        records.append({"Name": name, "CID": None})

results_df = pd.DataFrame(records)
results_df.to_csv("metabolite_cids.csv", index=False)
print(f"done: {len(results_df)}")