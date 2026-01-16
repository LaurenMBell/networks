import pandas as pd 

df1 = pd.read_csv("pls_DSS_CTRL_filtered_l.csv")
df2 = pd.read_csv("pls_DSS_CTRL_filtered.csv")

#df2.columns = df2.columns.str.replace('p-values', 'p_value')
#df2.columns = df2.columns.str.replace('p-value', 'p_value')

df2 = df2[df1.columns]

df1.convert_dtypes()
df2.convert_dtypes()

#df1 = df1.round(7)
#df2 = df2.round(7)

#df2.to_csv("lauren_fdr_plus_pval.csv")

count = 0

for line in range(len(df1)):
	if not df1.iloc[line].equals(df2.iloc[line]):
		count += 1
		print(f"MISMATCH AT LINE {line}")
		print(f"DF1:{df1.iloc[line]}")
		print(f"DF2:{df2.iloc[line]}\n") 

print(f"done, {count} mismatches between the two!!")

