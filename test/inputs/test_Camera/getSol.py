#%%
import numpy as np
import cv2
from scipy.io import loadmat

#%%
ncam = 5

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
def writeSol(file, pt_list, mode='w'):
    digit = "{:.8e}"
    with open(file, mode) as f:
        for i in range(ncam-1):
            pt = pt_list[i]
            size = pt.shape[0]
            for i in range(size):
                f.write(digit.format(pt[i,0]) + ',')
        
        pt = pt_list[ncam-1]
        for i in range(size-1):
            f.write(digit.format(pt[i,0]) + ',')
        f.write(digit.format(pt[size-1,0]) + '\n')

#%%
fileSol = '../../solutions/test_Camera/test_function_2.txt'

# get solution for projection
pt_world = np.array([1,2,3]).reshape(3,1)

# world to undistorted image
pt_camlist = []
for i in range(ncam):
    # project to camera space
    pt_cam = rotMatList[i] @ pt_world + transVecList[i]
    pt_cam = pt_cam / pt_cam[2]
    pt_camlist.append(pt_cam)
writeSol(fileSol, pt_camlist, 'w')


# undistorted image to distorted image  
pt_distlist = []
for i in range(ncam):
    pt_dist, _ = cv2.projectPoints(pt_world.reshape(1,1,3).astype(np.float32), rotVecList[i], transVecList[i], camMatList[i], distCoeffList[i])
    pt_distlist.append(pt_dist.reshape(2,1))
writeSol(fileSol, pt_distlist, 'a')

writeSol(fileSol, pt_distlist, 'a')

# %%
fileSol = '../../solutions/test_Camera/test_function_3.txt'

# get solution for calculate line of sight
pt_dist = np.array([461,441])

# distorted image to undistorted image
pt_undistlist = []
for i in range(ncam):
    pt_undist = cv2.undistortPointsIter(pt_dist.reshape(1,1,2).astype(np.float32), camMatList[i], distCoeffList[i], R=np.eye(3), P=camMatList[i], criteria=(cv2.TERM_CRITERIA_COUNT | cv2.TERM_CRITERIA_EPS, 50, 1e-5))
    
    pt_undist = pt_undist.reshape(2,1)
    pt_undist[0,0] = (pt_undist[0,0] - camMatList[i][0,2]) / camMatList[i][0,0]
    pt_undist[1,0] = (pt_undist[1,0] - camMatList[i][1,2]) / camMatList[i][1,1]
    
    pt_undistlist.append(pt_undist)
writeSol(fileSol, pt_undistlist, 'w')

# get line of sight 
pt_line_list = []
for i in range(ncam):
    pt = np.ones((3,1))
    pt[0,0] = pt_undistlist[i][0,0]
    pt[1,0] = pt_undistlist[i][1,0]
    
    pt_world = rotMatInvList[i] @ pt + transVecInvList[i]
    unit_vec = pt_world - transVecInvList[i]
    unit_vec = unit_vec / np.linalg.norm(unit_vec)
    
    res = np.zeros((6,1))
    res[:3,0] = transVecInvList[i].reshape(3)
    res[3:,0] = unit_vec.reshape(3)
    
    pt_line_list.append(res)
writeSol(fileSol, pt_line_list, 'a')

writeSol(fileSol, pt_line_list, 'a')



#%%
# Polynomial model
ncam = 2

err_list = []
ref_plane_list = []
n_coeff_list = []
u_coeff_list = []
v_coeff_list = []

for i in range(ncam):
    file = 'cam' + str(i+1) + '_poly' + '.txt'
    with open(file, 'r') as f:
        lines = f.readlines()[2:]
        if 'None' in lines[1] or 'none' in lines[1]:
            err_list.append(None)
        else:
            err_list.append(float(lines[1]))

        ref_plane = []
        content = lines[3].split(',')
        ref_plane.append(content[0])
        ref_plane.append(np.double(content[1]))
        ref_plane.append(np.double(content[2]))
        ref_plane_list.append(ref_plane)
        
        n_coeff = int(lines[5])
        n_coeff_list.append(n_coeff)
        
        u_coeff = np.zeros((n_coeff,4))
        for i in range(n_coeff):
            u_coeff[i,:] = np.array(lines[7+i].split(',')).astype(np.double)
        u_coeff_list.append(u_coeff)
        
        v_coeff = np.zeros((n_coeff,4))
        for i in range(n_coeff):
            v_coeff[i,:] = np.array(lines[7+n_coeff+1+i].split(',')).astype(np.double)
        v_coeff_list.append(v_coeff)


# %%
fileSol = '../../solutions/test_Camera/test_function_5.txt'

pt_world = np.array([1,0,3]).reshape(3,1)
# pt_world = np.array([1,0,3]).reshape(3,1)

pt_img_list = []
for camid in range(ncam):
    u = 0
    v = 0
    for i in range(n_coeff_list[camid]):
        u_coeff = u_coeff_list[camid][i, :]
        v_coeff = v_coeff_list[camid][i, :]
        u += u_coeff[0] * np.power(pt_world[0], int(u_coeff[1])) * np.power(pt_world[1], int(u_coeff[2])) * np.power(pt_world[2], int(u_coeff[3]))
        v += v_coeff[0] * np.power(pt_world[0], int(v_coeff[1])) * np.power(pt_world[1], int(v_coeff[2])) * np.power(pt_world[2], int(v_coeff[3]))
    pt_img = np.array([u,v]).reshape(2,1)
    pt_img_list.append(pt_img)
writeSol(fileSol, pt_img_list, 'w')

writeSol(fileSol, pt_img_list, 'a')
        
        


# %%
# load polynomial camera model from .mat
file = loadmat('polynomial.mat')
data = file['para_set']

cam3 = data[0,0]
cam3 = cam3[:int(cam3.shape[0]/2),:]

cam4 = data[0,1]
cam4 = cam4[:int(cam4.shape[0]/2),:]

# save to file 
def writePolyFile(path, cam):
    with open(path, 'w') as f:
        f.write('# Camera Model: (PINHOLE/POLYNOMIAL)\nPOLYNOMIAL\n')
        f.write('# Calibration Error: \nNone\n')
        f.write('# Reference Plane: (REF_X/REF_Y/REF_Z,coordinate,coordinate)\n')
        f.write('REF_Z,' + str(0) + ',' + str(5) + '\n')
        f.write('# Number of Coefficients: \n' + str(int(cam.shape[0]/2)) + '\n')
        f.write('# U_Coeff,X_Power,Y_Power,Z_Power\n')
        for i in range(int(cam.shape[0]/2)):
            f.write(str(cam[i,0]) + ',' + str(cam[i,1]) + ',' + str(cam[i,2]) + ',' + str(cam[i,3]) + '\n')
        f.write('# V_Coeff,X_Power,Y_Power,Z_Power\n')
        for i in range(int(cam.shape[0]/2), cam.shape[0]):
            f.write(str(cam[i,0]) + ',' + str(cam[i,1]) + ',' + str(cam[i,2]) + ',' + str(cam[i,3]) + '\n')

writePolyFile('cam3_poly.txt', cam3)
writePolyFile('cam4_poly.txt', cam4)


# %%
