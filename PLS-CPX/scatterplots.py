import pandas as pd
import matplotlib.pyplot as plt
import itertools
import numpy as np
from scipy.stats import spearmanr

def pairwise_scatterplot(gene, met, df):
    
    x = df.loc[gene].values.astype(float)
    y = df.loc[met].values.astype(float)

    plt.scatter(x, y)
    plt.xlabel(f"Gene: {gene}")
    plt.ylabel(f"Metabolite: {met}")

    z = np.polyfit(x, y, 1)
    pfit = np.poly1d(z)
    plt.plot(x, pfit(x), "r--")
    r, p = spearmanr(x, y)
    plt.text(0.05, 0.95, f"r = {r:.3f}\np = {p:.3e}",
             transform=plt.gca().transAxes,
            fontsize=12,
            verticalalignment='top')

    plt.show()
    return plt



lps = pd.read_csv("merged_lps.csv")
vecpac = pd.read_csv("merged_vecpac.csv")
dss = pd.read_csv("merged_dss.csv")

#rename columns to distinguish samples by model
dss_renamed = dss.rename(columns={col: f"{col}_DSS" if col != 'ID' else col for col in dss.columns})
vecpac_renamed = vecpac.rename(columns={col: f"{col}_VECPAC" if col != 'ID' else col for col in vecpac.columns})
lps_renamed = lps.rename(columns={col: f"{col}_LPS" if col != 'ID' else col for col in lps.columns})

#drop extra ID columns
vecpac_renamed.drop("ID", axis=1, inplace=True)
lps_renamed.drop("ID", axis=1, inplace=True)

#concatenate all models
all_data = pd.concat([dss_renamed, vecpac_renamed, lps_renamed], axis=1)




pairwise_scatterplot("ENSMUSG00000105302-P", "Trimethylsilyl 3-((trimethylsilyl)thio)propanoate", all_data)
#Expected: r=-0.044117647	p=0.871122005	
#ACTUAL: r=0.109 p=6.879e-01
pairwise_scatterplot("ENSMUSG00000027680-P", "Docosyl octyl ether", all_data)
#Expected: r=-0.356435859	p=0.210967224
#ACTUAL:  r=-0.218 p=4.545e-01
pairwise_scatterplot("ENSMUSG00000035293-P", "3,5-di-tert-Butyl-4-hydroxyacetophenone", all_data)
#expected: r=0.299862732	p=0.319547812
#actual: r=0.008 p=9.785e-01
pairwise_scatterplot("ENSMUSG00000019933-P", "Methylmalonic monoamide", all_data)
#expected: r=0.094609518	p=0.747669032
#actual: r = 0.231 p=4.264e-01
pairwise_scatterplot("ENSMUSG00000037110-P", "L-Proline", all_data)
#expected: r=0.28508404	p=0.284509067
#actual: r=0.254 3.423e-01