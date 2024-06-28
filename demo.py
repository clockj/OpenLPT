#%%
import numpy as np 
import matplotlib.pyplot as plt

import pyOpenLPT

# %%
out = pyOpenLPT.math.Matrix(2,2,3)
out[0] = 20

pyOpenLPT.math.matrix_to_numpy(out)
# %%
