#%%
import numpy as np 
import matplotlib.pyplot as plt

import pyOpenLPT
pyOpenLPT.PythonStreamRedirector()

#%%
# initialize a matrix
mtx_1 = pyOpenLPT.math.Matrix(2,2,0)
for i in range(4):
    mtx_1[i] = i
mtx_1.print(3)

# convert matrix to numpy array
mtx_2 = pyOpenLPT.math.matrix_to_numpy(mtx_1)
print(mtx_2)

# convert numpy array to matrix
mtx_3 = np.array([[1,2],[3,4]])
mtx_4 = pyOpenLPT.math.numpy_to_matrix(mtx_3)
mtx_4.print(3)


# initialize an image
img_1 = pyOpenLPT.math.Image(3,3,10)
img_1[1,1] = 1
img_1.print(3)

img_1_np = pyOpenLPT.math.matrix_to_numpy(img_1)
print(img_1_np)


# initialize Pt3D
pt3d = pyOpenLPT.math.Pt3D(1,2,3)
pt3d.print(3)
pt3d_np = pyOpenLPT.math.matrix_to_numpy(pt3d)
print(pt3d_np)

# initialize line3D
pt3d_1 = pyOpenLPT.math.Pt3D(1,2,3)
unit_vector = pyOpenLPT.math.Pt3D(1,1,1)
line3d = pyOpenLPT.math.Line3D()
line3d.pt = pt3d_1
line3d.unit_vector = unit_vector
line3d.pt.print(3)
line3d.unit_vector.print(3)


# matrix operations
mtx_5_np = np.array([[1,2],[3,4]])
mtx_5 = pyOpenLPT.math.numpy_to_matrix(mtx_5_np)

mtx_6_np = np.array([[5,6],[7,8]])
mtx_6 = pyOpenLPT.math.numpy_to_matrix(mtx_6_np)

mtx_sum = mtx_5 + mtx_6
print(mtx_sum==pyOpenLPT.math.numpy_to_matrix(mtx_5_np+mtx_6_np))

mtx_mul = mtx_5 * mtx_6
print(mtx_mul==pyOpenLPT.math.numpy_to_matrix(mtx_5_np@mtx_6_np))


mtx_sum += mtx_mul
mtx_sum.print(3)

mtx_sum += 1
mtx_sum.print(3)

mtx_sum -= 1
mtx_sum.print(3)

mtx_sum = mtx_5 + 1
mtx_sum.print(3)

mtx_sum *= 2
mtx_sum.print(3)

mtx_sum /= 2
mtx_sum.print(3)


pt3d_np = np.array([1,2,3]).reshape(3,1)
pt3d = pyOpenLPT.math.Pt3D(1,2,3)

rot_np = np.array([[1,0,0],[0,0,-1],[0,1,0]])
rot = pyOpenLPT.math.numpy_to_matrix(rot_np)

pt3d_rot = rot * pt3d

print(pt3d_rot==pyOpenLPT.math.numpy_to_matrix(rot_np@pt3d_np))

pt3d.T.print(3)
pt3d.transpose().print(3)
print(pt3d.norm())

# %%
# ImgIO 
imgIO = pyOpenLPT.math.ImageIO()
imgIO.loadImgPath('../../test/inputs/test_ImageIO/', 'test_function_1.txt')
img = imgIO.loadImg(0)

img = pyOpenLPT.math.matrix_to_numpy(img)
plt.imshow(img)

# %%
# Camera 
cam = pyOpenLPT.math.Camera('../../test/inputs/test_Camera/cam1.txt')
param = cam._pinhole_param

param.r_mtx.print(3)

# %%
