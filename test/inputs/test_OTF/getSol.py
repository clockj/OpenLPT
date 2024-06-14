#%%
import numpy as np


# %%
n_cam = 4
nx, ny, nz = 10, 10, 10
n_grid = nx * ny * nz

xmin, xmax, ymin, ymax, zmin, zmax = 0, 1, 0, 1, 0, 1

a = np.ones((n_cam, n_grid)) * 125
b = np.ones((n_cam, n_grid)) * 1.5
c = np.ones((n_cam, n_grid)) * 1.5
alpha = np.ones((n_cam, n_grid)) * 0

with open('otf.txt','w') as f:
    f.write('# Size: (n_cam,nx,ny,nz,n_grid)\n')
    f.write(f'{n_cam},{nx},{ny},{nz},{n_grid}\n')
    f.write('# Boundary: (xmin,xmax,ymin,ymax,zmin,zmax)\n')
    f.write(f'{xmin},{xmax},{ymin},{ymax},{zmin},{zmax}\n')
    f.write('# a: (n_cam*n_grid)\n')
    f.write('\n'.join([','.join(map(str, a[i])) for i in range(n_cam)]))
    f.write('\n')
    f.write('# b: (n_cam*n_grid)\n')
    f.write('\n'.join([','.join(map(str, b[i])) for i in range(n_cam)]))
    f.write('\n')
    f.write('# c: (n_cam*n_grid)\n')
    f.write('\n'.join([','.join(map(str, c[i])) for i in range(n_cam)]))
    f.write('\n')
    f.write('# alpha: (n_cam*n_grid)\n')
    f.write('\n'.join([','.join(map(str, alpha[i])) for i in range(n_cam)]))
    f.write('\n')
    

# %%
