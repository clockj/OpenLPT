#%%
import numpy as np 
# %%
R0 = np.array([1e-6,10e-6,100e-6,1e-3])
omega = 1/np.sqrt(1e3*R0*R0) * np.sqrt( 3*1.25*287*298 - 2*0.073/R0 )
omega2 = np.sqrt(3/(1e3*R0*R0) * ((1-0.0313)*101325 + 4*0.073/(3*R0)))
# %%
