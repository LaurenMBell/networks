import pandas as pd
import numpy as np

# Load both CSV files
mine = pd.read_csv("PLS-CPX_FDR.csv")
labmate = pd.read_csv("final_edges_with_fdr.csv")  # Replace with actual filename

# Drop rows with NaN in labmate's file
labmate = labmate.dropna()

print(f"My edges: {len(mine)}")
print(f"Labmate's edges after dropping NaN: {len(labmate)}")

# Rearrange your columns to match labmate's order
# Labmate's order: Gene, Metabolite, VECPAC r, VECPAC p, VECPAC n, DSS r, DSS p, DSS n, 
#                  LPS r, LPS p, LPS n, Pooled r, Pooled p, Pooled n, Consistent, Sign, Pooled FDR

mine_reordered = mine[[
    'gene',                    # Gene
    'metabolite',              # Metabolite
    'VECPAC r',               # VECPAC r
    'VECPAC p-values',        # VECPAC p
    'VECPAC n',               # VECPAC n
    'DSS r',                  # DSS r
    'DSS p-values',           # DSS p
    'DSS n',                  # DSS n
    'LPS r',                  # LPS r
    'LPS p-values',           # LPS p
    'LPS n',                  # LPS n
    'Pooled r',               # Pooled r
    'Pooled p-values',        # Pooled p
    'Pooled n',               # Pooled n
    'Consistent',             # Consistent
    'Sign',                   # Sign
    'Pooled FDR'              # Pooled FDR
]].copy()

# Rename columns to exactly match labmate's
mine_reordered.columns = [
    'Gene', 'Metabolite', 'VECPAC r', 'VECPAC p', 'VECPAC n',
    'DSS r', 'DSS p', 'DSS n', 'LPS r', 'LPS p', 'LPS n',
    'Pooled r', 'Pooled p', 'Pooled n', 'Consistent', 'Sign', 'Pooled FDR'
]

# Sort both dataframes by Gene and Metabolite (alphabetically)
mine_sorted = mine_reordered.sort_values(by=['Gene', 'Metabolite']).reset_index(drop=True)
labmate_sorted = labmate.sort_values(by=['Gene', 'Metabolite']).reset_index(drop=True)

# Save the rearranged and sorted files
mine_sorted.to_csv("PLS-CPX_FDR_lauren.csv", index=False)
labmate_sorted.to_csv("PLS-CPX_FDR_ananya.csv", index=False)

print(f"Saved sorted files:")
print(f"  My_FDR_Sorted.csv - {len(mine_sorted)} edges")
print(f"  Labmate_FDR_Sorted.csv - {len(labmate_sorted)} edges")

# Optional: Check if the gene-metabolite pairs match
my_pairs = set(zip(mine_sorted['Gene'], mine_sorted['Metabolite']))
labmate_pairs = set(zip(labmate_sorted['Gene'], labmate_sorted['Metabolite']))

common_pairs = my_pairs & labmate_pairs
only_mine = my_pairs - labmate_pairs
only_labmate = labmate_pairs - my_pairs

print(f"\nPair comparison:")
print(f"  Common pairs: {len(common_pairs)}")
print(f"  Only in mine: {len(only_mine)}")
print(f"  Only in labmate's: {len(only_labmate)}")