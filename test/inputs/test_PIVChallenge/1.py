#%%
import numpy as np
import matplotlib.pyplot as plt
import cv2
import os
# %%
calib_data = np.loadtxt('D_cal/D_cal.csv', delimiter=',', skiprows=1)
width = 4512
height = 800
imgSize = (width, height)
n_cam = int((calib_data.shape[1]-3)/2)

# %%
plt.figure()

cam_id = 3
plt.plot(calib_data[:, 3+cam_id*2], calib_data[:, 4+cam_id*2], 'r.')

#%%
# pinhole model for plate calibration
def calibPinhole (pt3d, pt2d, plane_direction='z'):
    axis = 2
    if plane_direction == 'y':
        axis = 1
    elif plane_direction == 'x':
        axis = 0
    
    # split points by plane
    plane_id = np.unique(pt3d[:, axis])
    
    # check if 0 is included in plane_id
    if 0 not in plane_id:
        print('Plane 0 is not included in the calibration data.')
        return -1, None, None, None, None
    n_plane = len(plane_id)
    
    # check the data type
    pt3d = pt3d.astype(np.float32)
    pt2d = pt2d.astype(np.float32)
    
    # get camMat based on plane_id=0
    judge = pt3d[:, axis] == 0
    pt3d_0 = pt3d[judge,:].reshape(-1,1,3)
    pt2d_0 = pt2d[judge,:].reshape(-1,1,2)
    _, camMat_0, _, _, _ = cv2.calibrateCamera(
        [pt3d_0], [pt2d_0], imgSize, None, None
    )
    
    pt3d_list = []
    pt2d_list = []
    for i in range(n_plane):
        judge = pt3d[:, axis] == plane_id[i]
        pt3d_list.append(pt3d[judge,:].astype(np.float32).reshape(-1,1,3))
        pt2d_list.append(pt2d[judge,:].astype(np.float32).reshape(-1,1,2))
        
    # calibrate camMat and distCoeff
    camCalibErr, camMat, distCoeff, _, _ = cv2.calibrateCamera(
        pt3d_list, pt2d_list, imgSize, camMat_0, None, flags=cv2.CALIB_USE_INTRINSIC_GUESS + cv2.CALIB_FIX_K1 + cv2.CALIB_FIX_K2 + cv2.CALIB_FIX_K3 + cv2.CALIB_ZERO_TANGENT_DIST
    )
    
    # calibrate rotation and translation
    optFlag = cv2.SOLVEPNP_ITERATIVE
    _, rotVec, transVec = cv2.solvePnP(
        pt3d.reshape(-1,1,3), pt2d.reshape(-1,1,2), camMat, distCoeff, flags=optFlag
    )
    
    # calculate error    
    pt2d_proj, _ = cv2.projectPoints(pt3d, rotVec, transVec, camMat, distCoeff)
    pt2d_proj = pt2d_proj.reshape(-1,2)
    poseCalibErr = cv2.norm(pt2d, pt2d_proj, cv2.NORM_L2) / pt2d.shape[0]
        
    return camCalibErr, poseCalibErr, camMat, distCoeff, rotVec, transVec 

# imgSize=(width, height)
def saveCamFile (filename, camCalibErr, poseCalibErr, imgSize, camMat, distCoeff, rotVec, transVec):
    with open(filename, 'w') as f:
        f.write('# Camera Model: (PINHOLE/POLYNOMIAL)\n' + str('PINHOLE') + '\n')
        f.write('# Camera Calibration Error: \n' + str(camCalibErr) + '\n')
        f.write('# Pose Calibration Error: \n' + str(poseCalibErr) + '\n')
        
        f.write('# Image Size: (n_row,n_col)\n')
        f.write(str(imgSize[1])+','+str(imgSize[0])+'\n') # OpenCV: imgSize=(width, height)
        
        f.write('# Camera Matrix: \n')
        f.write(str(camMat[0,0])+','+str(camMat[0,1])+','+str(camMat[0,2])+'\n')
        f.write(str(camMat[1,0])+','+str(camMat[1,1])+','+str(camMat[1,2])+'\n')
        f.write(str(camMat[2,0])+','+str(camMat[2,1])+','+str(camMat[2,2])+'\n')
        f.write('# Distortion Coefficients: \n')
        f.write(str(distCoeff[0,0])+','+str(distCoeff[0,1])+','+str(distCoeff[0,2])+','+str(distCoeff[0,3])+','+str(distCoeff[0,4])+'\n')
        f.write('# Rotation Vector: \n')
        f.write(str(rotVec[0,0])+','+str(rotVec[1,0])+','+str(rotVec[2,0])+'\n')
        f.write('# Rotation Matrix: \n')
        rotMat = cv2.Rodrigues(rotVec)[0]
        f.write(str(rotMat[0,0])+','+str(rotMat[0,1])+','+str(rotMat[0,2])+'\n')
        f.write(str(rotMat[1,0])+','+str(rotMat[1,1])+','+str(rotMat[1,2])+'\n')
        f.write(str(rotMat[2,0])+','+str(rotMat[2,1])+','+str(rotMat[2,2])+'\n')
        f.write('# Inverse of Rotation Matrix: \n')
        rotMatInv = np.linalg.inv(rotMat)
        f.write(str(rotMatInv[0,0])+','+str(rotMatInv[0,1])+','+str(rotMatInv[0,2])+'\n')
        f.write(str(rotMatInv[1,0])+','+str(rotMatInv[1,1])+','+str(rotMatInv[1,2])+'\n')
        f.write(str(rotMatInv[2,0])+','+str(rotMatInv[2,1])+','+str(rotMatInv[2,2])+'\n')
        f.write('# Translation Vector: \n')
        f.write(str(transVec[0,0])+','+str(transVec[1,0])+','+str(transVec[2,0])+'\n')
        f.write('# Inverse of Translation Vector: \n')
        transVecInv = -np.matmul(rotMatInv, transVec)
        f.write(str(transVecInv[0,0])+','+str(transVecInv[1,0])+','+str(transVecInv[2,0])+'\n')

def saveImgFile (filename, frame_start, frame_end, format):
    with open(filename, 'w') as f:
        for i in range(frame_start, frame_end+1):
            f.write(format.format(i) + '\n')

# %%
# Calibration
camCalibErrList = []
poseCalibErrList = []
camMatList = []
distCoeffList = []
rotVecList = []
rotMatList = []
transVecList = []

# Pinhole model 
distCoeffs = np.zeros(5)
pt3d = calib_data[:, 0:3]
for i in range(n_cam):
    pt2d = calib_data[:, 3+i*2:5+i*2]
    camCalibErr, poseCalibErr, camMat, distCoeff, rotVec, transVec = calibPinhole (pt3d, pt2d, plane_direction='z')
    rotMat = cv2.Rodrigues(rotVec)[0]
    
    camCalibErrList.append(camCalibErr)
    poseCalibErrList.append(poseCalibErr)
    camMatList.append(camMat)
    distCoeffList.append(distCoeff)
    rotVecList.append(rotVec)
    rotMatList.append(rotMat)
    transVecList.append(transVec)
    
    saveCamFile('camFile/cam'+str(i+1)+'.txt', camCalibErr, poseCalibErr, imgSize, camMat, distCoeff, rotVec, transVec)

#%%
img = cv2.imread('D_Img/D_cam3_0001.tif',cv2.IMREAD_ANYDEPTH)
plt.imshow(img.astype(np.float32)/(2**8-1), cmap='gray', vmax=10)

# %%
# save image file 
path_curr = '..\\test\\inputs\\test_PIVChallenge\\'

for i in range(n_cam):
    path_img = os.path.join(path_curr, 'D_Img', 'D_cam'+str(i)+'_'+'{:04d}.tif')
    
    saveImgFile('imgFile/cam'+str(i+1)+'ImageNames.txt', 0, 50, path_img)

# %%
volume_size = np.array([204.8, 25.6, 17.6])
limit = volume_size/2

print(limit)
# %%
