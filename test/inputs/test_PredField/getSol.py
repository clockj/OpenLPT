#%%
import numpy as np


#%%
np.random.seed(1234)

npts = int(1.5e4)
pts = np.random.rand(npts, 3) * 40 - 20 # 3D points in [-20, 20]

# suppose nx,ny,nz = 51,51,51 
# gird size = 40/50 = 0.8
# search radius = 0.8 / 2 = 0.4
n_xyz = [51, 51, 51]
n_grid = np.prod(n_xyz)
grid = np.zeros((n_grid, 3))

grid_x = np.linspace(-20, 20, n_xyz[0])
grid_y = np.linspace(-20, 20, n_xyz[1])
grid_z = np.linspace(-20, 20, n_xyz[2])

id = 0
for i in range(n_xyz[0]):
    for j in range(n_xyz[1]):
        for k in range(n_xyz[2]):
            grid[id] = [grid_x[i], grid_y[j], grid_z[k]]
            id += 1

def mapGridID (x_id, y_id, z_id, n_xyz):
    return x_id * n_xyz[1] * n_xyz[2] + y_id * n_xyz[2] + z_id
    
#%%
# const velocity
vel = 0.1

# calculate new position
pts_new = pts + vel

# calculate disp field
disp_field = vel * np.ones((n_grid, 3))

np.savetxt('pts_prev_1.csv', pts, fmt='%.8f', delimiter=',')
np.savetxt('pts_curr_1.csv', pts_new, fmt='%.8f', delimiter=',')
np.savetxt('../../solutions/test_PredField/disp_field_1.csv', disp_field, fmt='%.8f', delimiter=',')

# %%
# ideal vortex: vorticity direction is along z-axis
gamma = 0.05

# calculate new position
radius = np.linalg.norm(pts[:,:2], axis=1)
unit_vec = pts[:,:2] / radius.reshape(-1,1)
perp_vec = np.array([-unit_vec[:,1], unit_vec[:,0], np.zeros(npts)]).T

vel = gamma / radius.reshape(-1,1) * perp_vec
pts_new = pts + vel

# calculate disp field
radius = np.linalg.norm(grid, axis=1)
unit_vec = grid[:,:2] / radius.reshape(-1,1)
perp_vec = np.array([-unit_vec[:,1], unit_vec[:,0], np.zeros(n_grid)]).T

vel = gamma / radius.reshape(-1,1) * perp_vec
vel = np.nan_to_num(vel, nan=0)
disp_field = vel

# save 
np.savetxt('pts_prev_2.csv', pts, fmt='%.8f', delimiter=',')
np.savetxt('pts_curr_2.csv', pts_new, fmt='%.8f', delimiter=',')
np.savetxt('../../solutions/test_PredField/disp_field_2.csv', disp_field, fmt='%.8f', delimiter=',')

# %%
