import qnorm
import pandas as pd
import matplotlib.pyplot as plt

def validation_1(filename):
	df = pd.read_csv(f"Lauren_Analysis_files/pls_{filename}_qnormed.csv", header=[0,1])
	summary  = df.describe()
	print(summary)
	

validation_1("DSS")
validation_1("LPS")
validation_1("VECPAC")

def validation_2():
	
	mtx = pd.DataFrame({'C1': {'A': 5, 'B': 2, 'C': 3, 'D': 4},
                   'C2': {'A': 4, 'B': 1, 'C': 4, 'D': 2},
                   'C3': {'A': 3, 'B': 4, 'C': 6, 'D': 8}})
	
	print("matrix before nomalization:\n")
	print(mtx)

	print("expected matrix after normalization:\n")
	print("C1        C2        C3\nA  5.666667  5.166667  2.000000\nB  2.000000  2.000000  3.000000\nC  3.000000  5.166667  4.666667\nD  4.666667  3.000000  5.666667")
	
	print("mactual matrix after normalization:\n")
	
	print(qnorm.quantile_normalize(mtx.astype(float), axis=1))
	
validation_2()
	
