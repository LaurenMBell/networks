import pandas as pd

m = input("enter 1 for single file cleaning, 2 for 2 file cleaning")
if m == 1:
    print("cleaning one file")
    
elif m==2:
    A = pd.read_csv("final_edges_with_fdr 2.csv").dropna()
    L = pd.read_csv("12-5_FDR_L.csv")

    A_clean = A[["Gene", "Metabolite", "Pooled r", "Pooled p", "Pooled FDR"]].copy()

    A_clean.columns = ["Gene", "Metabolite", "r", "p", "fdr"]


    A_clean["Gene"] = A_clean["Gene"].astype(str).str.replace("-P", "", regex=False)

    L_clean = L[["gene", "metabolite", "Pooled r", "Pooled p-values", "Pooled FDR"]].copy()

    L_clean.columns = ["Gene", "Metabolite", "r", "p", "fdr"]

    L_clean["Gene"] = L_clean["Gene"].astype(str).str.replace("-P", "", regex=False)

    A_clean_sorted = A_clean.sort_values(["Gene", "Metabolite"]).reset_index(drop=True)
    L_clean_sorted = L_clean.sort_values(["Gene", "Metabolite"]).reset_index(drop=True)

    A_clean_sorted.to_csv("12-5_FDR_A_reformatted.csv", index=False)
    L_clean_sorted.to_csv("12-5_FDR_L_reformatted.csv", index=False)
