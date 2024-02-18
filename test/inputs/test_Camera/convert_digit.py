#%%
import numpy as np
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
#%%
digit = "{:.8e}"
for i in range(ncam):
    file = '../../solutions/test_Camera/'+'cam'+str(i+1)+'.txt'
    camMat = camMatList[i]
    distCoeff = distCoeffList[i]
    rotVec = rotVecList[i]
    rotMat = rotMatList[i]
    rotMatInv = rotMatInvList[i]
    transVec = transVecList[i]
    transVecInv = transVecInvList[i]
    
    with open(file, 'w') as f:
        f.write('# Camera Model: (PINHOLE/POLYNOMIAL)\n' + str('PINHOLE') + '\n')
        f.write('# Camera Calibration Error: \n' + str(None) + '\n')
        f.write('# Pose Calibration Error: \n' + str(None) + '\n')
        
        f.write('# Camera Matrix: \n')
        for i in range(3):
            for j in range(2):
                f.write(digit.format(camMat[i,j]) + ',')
            f.write(digit.format(camMat[i,2]) + '\n')
            
        f.write('# Distortion Coefficients: \n')
        for i in range(4):
            f.write(digit.format(distCoeff[0,i]) + ',')
        f.write(digit.format(distCoeff[0,4]) + '\n')
        
        f.write('# Rotation Vector: \n')
        for i in range(2):
            f.write(digit.format(rotVec[i,0]) + ',')
        f.write(digit.format(rotVec[2,0]) + '\n')
        
        f.write('# Rotation Matrix: \n')
        for i in range(3):
            for j in range(2):
                f.write(digit.format(rotMat[i,j]) + ',')
            f.write(digit.format(rotMat[i,2]) + '\n')
        
        f.write('# Inverse of Rotation Matrix: \n')
        for i in range(3):
            for j in range(2):
                f.write(digit.format(rotMatInv[i,j]) + ',')
            f.write(digit.format(rotMatInv[i,2]) + '\n')
        
        f.write('# Translation Vector: \n')
        for i in range(2):
            f.write(digit.format(transVec[i,0]) + ',')
        f.write(digit.format(transVec[2,0]) + '\n')
        
        f.write('# Inverse of Translation Vector: \n')
        for i in range(2):
            f.write(digit.format(transVecInv[i,0]) + ',')
        f.write(digit.format(transVecInv[2,0]) + '\n')
# %%
