#%%
import numpy as np
import cv2
from scipy.io import loadmat
import matplotlib.pyplot as plt
import getTiffImg
#%%
ncam = 4

camcalibErrList = []
posecalibErrList = []
imgSizeList = []
camMatList = []
distCoeffList = []
rotVecList = []
rotMatList = []
rotMatInvList = []
transVecList = []
transVecInvList = []

# Load camera parameters
for i in range(ncam):
    file = 'cam' + str(i+1) + '.txt'
    with open(file, 'r') as f:
        line_id = 0
        
        lines = f.readlines()[2:]
        
        line_id += 1
        if 'None' in lines[line_id] or 'none' in lines[line_id]:
            camcalibErrList.append(None)
        else:
            camcalibErrList.append(float(lines[line_id]))
        line_id += 2
        if 'None' in lines[line_id] or 'none' in lines[line_id]:
            posecalibErrList.append(None)
        else:
            posecalibErrList.append(float(lines[line_id]))
        
        line_id += 2
        imgSize = np.array(lines[line_id].split(',')).astype(np.int32)
        imgSizeList.append(imgSize)
        
        line_id += 2
        camMat = np.zeros((3,3))
        camMat[0,:] = np.array(lines[line_id].split(',')).astype(np.double)
        camMat[1,:] = np.array(lines[line_id+1].split(',')).astype(np.double)
        camMat[2,:] = np.array(lines[line_id+2].split(',')).astype(np.double)
        camMatList.append(camMat)
        
        line_id += 4
        distCoeff = np.array([lines[line_id].split(',')]).astype(np.double)
        distCoeffList.append(distCoeff)
        
        line_id += 2
        rotVec = np.zeros((3,1))
        rotVec[:,0] = np.array(lines[line_id].split(',')).astype(np.double)
        rotVecList.append(rotVec)
        
        line_id += 2
        rotMat = np.zeros((3,3))
        rotMat[0,:] = np.array(lines[line_id].split(',')).astype(np.double)
        rotMat[1,:] = np.array(lines[line_id+1].split(',')).astype(np.double)
        rotMat[2,:] = np.array(lines[line_id+2].split(',')).astype(np.double)
        rotMatList.append(rotMat)
        
        # rotMatInv
        line_id += 4
        rotMatInv = np.zeros((3,3))
        rotMatInv[0,:] = np.array(lines[line_id].split(',')).astype(np.double)
        rotMatInv[1,:] = np.array(lines[line_id+1].split(',')).astype(np.double)
        rotMatInv[2,:] = np.array(lines[line_id+2].split(',')).astype(np.double)
        rotMatInvList.append(rotMatInv)
        
        line_id += 4
        transVec = np.zeros((3,1))
        transVec[:,0] = np.array(lines[line_id].split(',')).astype(np.double)
        transVecList.append(transVec)
        
        # transVecInv
        line_id += 2
        transVecInv = np.zeros((3,1))
        transVecInv[:,0] = np.array(lines[line_id].split(',')).astype(np.double)
        transVecInvList.append(transVecInv)

# %%
# randomly generate 3D points
# npts = int(1.5e4)
npts = int(4e4)

np.random.seed(1234)
pt3d_list = np.random.rand(npts, 3) * 40 - 20
is_select = np.ones(npts, dtype=bool)

for i in range(ncam):
    # project 3D points to 2D
    pt2d_list = cv2.projectPoints(pt3d_list.reshape(npts, 1, 3), rotVecList[i], transVecList[i], camMatList[i], distCoeffList[i])[0].reshape(npts, 2)

    # remove points outside the image
    is_select = np.logical_and(is_select, np.logical_and(pt2d_list[:,0] > 0, pt2d_list[:,0] < imgSizeList[i][1]))
    is_select = np.logical_and(is_select, np.logical_and(pt2d_list[:,1] > 0, pt2d_list[:,1] < imgSizeList[i][0]))

pt3d_list = pt3d_list[is_select]

pt2d_list_all = []
for i in range(ncam):
    pt2d_list = cv2.projectPoints(pt3d_list.reshape(-1, 1, 3), rotVecList[i], transVecList[i], camMatList[i], distCoeffList[i])[0].reshape(-1, 2)
    pt2d_list_all.append(pt2d_list)
pt2d_list_all = np.array(pt2d_list_all)


# save to csv file 
np.savetxt('../../solutions/test_StereoMatch/pt3d_list.csv', pt3d_list, delimiter=',', fmt='%.8f')

for i in range(ncam):
    np.savetxt('../../inputs/test_StereoMatch/pt2d_list_cam' + str(i+1) + '.csv', pt2d_list_all[i], delimiter=',', fmt='%.8f')
    
# %%
# generate tiff images

# randomly generate 3D points
npts = int(1.5e4)
# npts = int(4e4)

np.random.seed(1234)
pt3d_list = np.random.rand(npts, 3) * 40 - 20
is_select = np.ones(npts, dtype=bool)
for i in range(ncam):
    # project 3D points to 2D
    pt2d_list = cv2.projectPoints(pt3d_list.reshape(npts, 1, 3), rotVecList[i], transVecList[i], camMatList[i], distCoeffList[i])[0].reshape(npts, 2)
    # remove points outside the image
    is_select = np.logical_and(is_select, np.logical_and(pt2d_list[:,0] > 0, pt2d_list[:,0] < imgSizeList[i][1]))
    is_select = np.logical_and(is_select, np.logical_and(pt2d_list[:,1] > 0, pt2d_list[:,1] < imgSizeList[i][0]))

pt3d_list = pt3d_list[is_select]
pt2d_list_all = []
for i in range(ncam):
    pt2d_list = cv2.projectPoints(pt3d_list.reshape(-1, 1, 3), rotVecList[i], transVecList[i], camMatList[i], distCoeffList[i])[0].reshape(-1, 2)
    pt2d_list_all.append(pt2d_list)
pt2d_list_all = np.array(pt2d_list_all)

np.savetxt('../../solutions/test_StereoMatch/pt3d_list_img.csv', pt3d_list, delimiter=',', fmt='%.8f')

#%%
# generate tiff images
img_list = []
for i in range(ncam):
    img = getTiffImg.getTiffImg(pt3d_list, rotVecList[i], transVecList[i], camMatList[i], distCoeffList[i], imgSizeList[i])
    img_list.append(img)

    # save to tiff files
    cv2.imwrite('../../inputs/test_StereoMatch/img_cam'+str(i+1)+'.tif', img_list[i])

# %%
# load results 
tr3d_list = np.loadtxt('../../results/test_StereoMatch/tr3d_img_deleteGhost.csv', delimiter=',', skiprows=1)

pt3d_list_res = tr3d_list[:,:3]

npts = pt3d_list.shape[0]

tor = 5e-3
mismatchID_list = []
for i in range(npts):
    error = np.linalg.norm(pt3d_list[i] - pt3d_list_res, axis=1)
    min_id = np.argmin(error)
    if error[min_id] > tor:
        mismatchID_list.append(i)
        
n_mismatch = len(mismatchID_list)
rate_mismatch = n_mismatch / npts * 100

print('Number of mismatched points:', n_mismatch)
print('Mismatch rate:', rate_mismatch, '%')

# %%
