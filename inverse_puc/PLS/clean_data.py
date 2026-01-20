#note for lauren: this isn't apart of the pipeline 

import pandas as pd 

pls = pd.read_csv("PLS_edges.csv")
pls = pls[(pls["FDR"] <= 0.05) & (pls[["VECPAC p-values", "DSS p-values", "LPS p-values"]].max(axis=1) <= 0.2)]

pls_cpx = pd.read_csv("PLS-CPX_edges.csv")
pls_cpx = pls_cpx[pls_cpx["Pooled FDR"] <= 0.1]


pls.sort_values("Metabolite 1")
pls_cpx.sort_values("cpx_gene")

pls.to_csv("PLS_edges_sorted.csv", index=False)
pls_cpx.to_csv("PLS-CPX_edges_sorted.csv", index=False)
