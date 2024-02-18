#%%
import numpy as np
import cv2

#%%
ncam = 5

camcalibErrList = []
posecalibErrList = []
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
        lines = f.readlines()[2:]
        
        if 'None' in lines[1] or 'none' in lines[1]:
            camcalibErrList.append(None)
        else:
            camcalibErrList.append(float(lines[1]))
        if 'None' in lines[3] or 'none' in lines[3]:
            posecalibErrList.append(None)
        else:
            posecalibErrList.append(float(lines[3]))
        
        camMat = np.zeros((3,3))
        camMat[0,:] = np.array(lines[5].split(',')).astype(np.double)
        camMat[1,:] = np.array(lines[6].split(',')).astype(np.double)
        camMat[2,:] = np.array(lines[7].split(',')).astype(np.double)
        camMatList.append(camMat)
        
        distCoeff = np.array([lines[9].split(',')]).astype(np.double)
        distCoeffList.append(distCoeff)
        
        rotVec = np.zeros((3,1))
        rotVec[:,0] = np.array(lines[11].split(',')).astype(np.double)
        rotVecList.append(rotVec)
        
        rotMat = np.zeros((3,3))
        rotMat[0,:] = np.array(lines[13].split(',')).astype(np.double)
        rotMat[1,:] = np.array(lines[14].split(',')).astype(np.double)
        rotMat[2,:] = np.array(lines[15].split(',')).astype(np.double)
        rotMatList.append(rotMat)
        # line 17,18,19 are rotMatInv
        rotMatInv = np.zeros((3,3))
        rotMatInv[0,:] = np.array(lines[17].split(',')).astype(np.double)
        rotMatInv[1,:] = np.array(lines[18].split(',')).astype(np.double)
        rotMatInv[2,:] = np.array(lines[19].split(',')).astype(np.double)
        rotMatInvList.append(rotMatInv)
        
        transVec = np.zeros((3,1))
        transVec[:,0] = np.array(lines[21].split(',')).astype(np.double)
        transVecList.append(transVec)
        # line 23 is transVecInv
        transVecInv = np.zeros((3,1))
        transVecInv[:,0] = np.array(lines[23].split(',')).astype(np.double)
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

# %%
