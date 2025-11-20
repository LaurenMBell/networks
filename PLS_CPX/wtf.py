import pandas as pd
table1 = pd.read_csv("PLS-CPX_FDR_matt.csv")
table2 = pd.read_csv("PLS-CPX_FDR_reordered.csv")

table1.loc['n1'] = table1.loc['n1'].strip("-P")

table1_gene = set(table1['n1'])
table2_gene = set(table2['n2'])

print(table1_gene.difference(table2_gene))