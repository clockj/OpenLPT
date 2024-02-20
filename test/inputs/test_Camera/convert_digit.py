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
#%%
digit = "{:.8e}"
for cam_id in range(ncam):
    file = '../../solutions/test_Camera/'+'cam'+str(cam_id+1)+'_poly'+'.txt'
    with open(file, 'w') as f:
        f.write('# Camera Model: (PINHOLE/POLYNOMIAL)\n' + str('POLYNOMIAL') + '\n')
        f.write('# Camera Calibration Error: \n' + str(None) + '\n')
        
        f.write('# Reference Plane: (REF_X/REF_Y/REF_Z,coordinate,coordinate)\n')
        f.write(ref_plane_list[cam_id][0] + ',' + str(int(ref_plane_list[cam_id][1])) + ',' + str(int(ref_plane_list[cam_id][2])) + '\n')
        
        f.write('# Number of Coefficients: \n' + str(n_coeff_list[cam_id]) + '\n')
        
        f.write('# U_Coeff,X_Power,Y_Power,Z_Power\n')
        for i in range(n_coeff_list[cam_id]):
            for j in range(3):
                f.write(digit.format(u_coeff_list[cam_id][i,j]) + ',')
            f.write(digit.format(u_coeff_list[cam_id][i,3]) + '\n')
        
        f.write('# V_Coeff,X_Power,Y_Power,Z_Power\n')
        for i in range(n_coeff_list[cam_id]):
            for j in range(3):
                f.write(digit.format(v_coeff_list[cam_id][i,j]) + ',')
            f.write(digit.format(v_coeff_list[cam_id][i,3]) + '\n')
    

# %%
