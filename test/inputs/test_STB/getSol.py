#%%
import numpy as np
import cv2
from scipy.io import loadmat
import matplotlib.pyplot as plt
import getTiffImg
import os

# %%
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
folder = 'camFile/'
for i in range(ncam):
    file = folder + 'cam' + str(i+1) + '.txt'
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
# load tracks 
tracks = loadmat('tracks_12k5_coarse_250frames.mat')['tracks']

print([np.min(tracks[:,3]),np.max(tracks[:,3])])
print([np.min(tracks[:,4]),np.max(tracks[:,4])])

# %%
# check the number of cameras that can be seen for each point 
pts = tracks[:,0:3]
npts = pts.shape[0]
npts_cam = np.zeros(npts, dtype=np.int32)

for i in range(ncam):
    pts_proj = cv2.projectPoints(pts.reshape(npts,1,3), rotVecList[i], transVecList[i], camMatList[i], distCoeffList[i])[0].reshape(npts, 2)

    judge = np.all((pts_proj[:,0] >= 1, pts_proj[:,0] < imgSizeList[i][1]-1, pts_proj[:,1] >= 1, pts_proj[:,1] < imgSizeList[i][0]-1), axis=0)
    
    npts_cam[judge] += 1
    
#%%

    
# %%
# generate folders
folder = 'imgFile/'

for i in range(ncam):
    dir = folder+'cam'+str(i+1)+'/'
    if not os.path.exists(dir):
        os.makedirs(dir)

# frame range 
frame_range = [0, 49]

# generate image file names
format_str = '../test/inputs/test_STB/imgFile/cam{:d}/img{:05d}.tif\n'
for i in range(ncam):
    file = folder + 'cam' + str(i+1) + 'ImageNames.txt'
    
    with open(file, 'w') as f:
        for j in range(frame_range[0], frame_range[1]+1):
            f.write(format_str.format(i+1, j))
        

# # generate tiff images
# for i in range(frame_range[0], frame_range[1]+1):    
#     for j in range(ncam):
#         img = getTiffImg.getTiffImg(tracks[tracks[:,3]==i+1,0:3], rotVecList[j], transVecList[j], camMatList[j], distCoeffList[j], imgSizeList[j])

#         file = folder + 'cam' + str(j+1) + '/img' + '{:05d}'.format(i) + '.tif'
#         cv2.imwrite(file, img)
# %%
