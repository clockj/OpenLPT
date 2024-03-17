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
# vortex: vorticity direction is along z-axis
def vortex(pts):
    gamma = 0.05
    r0 = 1
    
    vel = np.zeros_like(pts)
    radius = np.linalg.norm(pts[:,:2], axis=1)
    is_linear = radius < r0 
    is_vortex = np.bitwise_not(is_linear)

    unit_vec = pts[:,:2] / radius.reshape(-1,1)
    unit_vec = np.nan_to_num(unit_vec, nan=0)
    perp_vec = np.array([-unit_vec[:,1], unit_vec[:,0], np.zeros_like(unit_vec[:,0])]).T
    
    vel[is_linear,:] = gamma * radius[is_linear].reshape(-1,1) * perp_vec[is_linear]

    vel[is_vortex,:] = gamma / radius[is_vortex].reshape(-1,1) * perp_vec[is_vortex]
    
    return vel

# calculate new position
vel = vortex(pts)
pts_new = pts + vel

# calculate disp field
disp_field = vortex(grid)

# save 
np.savetxt('pts_prev_2.csv', pts, fmt='%.8f', delimiter=',')
np.savetxt('pts_curr_2.csv', pts_new, fmt='%.8f', delimiter=',')
np.savetxt('../../solutions/test_PredField/disp_field_2.csv', disp_field, fmt='%.8f', delimiter=',')

# %%
