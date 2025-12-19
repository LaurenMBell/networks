from scipy import stats
import numpy as np
tnan = [1, 2, 3, 4, 5, np.nan]
unan = [0, 0, 0, 0, 0, np.nan]
r, p = stats.spearmanr(tnan, unan, nan_policy='omit')

print(r)
print(p)
