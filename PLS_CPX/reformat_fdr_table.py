import pandas as pd
import numpy as np

table1 = pd.read_csv("PLS-CPX_FDR.csv")
table2 = pd.read_csv("fdr_table.csv") 
table2 = table2.dropna()

"""""
table1_reordered = table1[[
    'gene',                    
    'metabolite',              
    #'VECPAC r',              
    #'VECPAC p-values',       
    #'VECPAC n',               
    #'DSS r',                 
    #'DSS p-values',          
    #'DSS n',                  
    #'LPS r',                  
    #'LPS p-values',           
    #'LPS n',                 
    'Pooled r',              
    'Pooled p-values',       
    #'Pooled n',               
    #'Consistent',             
    #'Sign',                   
    'Pooled FDR'              
]].copy()
"""

table1.columns = [
#table1_reordered.columns = [
    'n1', 'n2', #'VECPAC r', 'VECPAC p', 'VECPAC n',
    #'DSS r', 'DSS p', 'DSS n', 'LPS r', 'LPS p', 'LPS n',
    'r', 'p', 'fdr'#'Pooled n', 'Consistent', 'Sign', 'Pooled FDR'
]

table1_sorted = table1.sort_values(by=['n1', 'n2']).reset_index(drop=True)
table2_sorted = table2.sort_values(by=['n1', 'n2']).reset_index(drop=True)

table1_sorted.to_csv("PLS-CPX_FDR.csv", index=False)
table2_sorted.to_csv("PLS-CPX_FDR_matt.csv", index=False)
