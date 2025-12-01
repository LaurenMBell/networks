import pandas as pd
import numpy as np

table2 = pd.read_csv("PLS-CPX_FDR.csv")   
table1 = pd.read_csv("fdr_table.csv").dropna()  


table2_reordered = table2[[
    'gene',
    'metabolite',
    'Pooled r',
    'Pooled p-values',
    'Pooled FDR'
]].copy()

table2_reordered.columns = ['n1', 'n2', 'r', 'p', 'fdr']
table2_reordered['n1'] = table2_reordered['n1'].str.replace('-P', '', regex=False)

table1_sorted = table1.sort_values(by=['n1', 'n2']).reset_index(drop=True)
table2_sorted = table2_reordered.sort_values(by=['n1', 'n2']).reset_index(drop=True)

table1_sorted.to_csv("PLS-CPX_FDR_matt.csv", index=False)
table2_sorted.to_csv("PLS-CPX_FDR_reordered.csv", index=False)
