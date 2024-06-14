#%%
import numpy as np
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

#%%
digit = "{:.8e}"
for i in range(ncam):
    file = '../../solutions/test_Camera/'+'cam'+str(i+1)+'.txt'
    imgSize = imgSizeList[i]
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
        
        f.write('# Image Size: (n_row,n_col)\n')
        f.write(str(imgSize[0]) + ',' + str(imgSize[1]) + '\n')
        
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
imgSize_list = []
ref_plane_list = []
n_coeff_list = []
u_coeff_list = []
v_coeff_list = []

for i in range(ncam):
    file = 'cam' + str(i+1) + '_poly' + '.txt'
    with open(file, 'r') as f:
        lines = f.readlines()[2:]
        
        line_id = 0
        line_id += 1
        if 'None' in lines[line_id] or 'none' in lines[line_id]:
            err_list.append(None)
        else:
            err_list.append(float(lines[line_id]))

        line_id += 2
        imgSize = np.array(lines[line_id].split(',')).astype(np.int32)
        imgSize_list.append(imgSize)
        
        line_id += 2
        ref_plane = []
        content = lines[line_id].split(',')
        ref_plane.append(content[0])
        ref_plane.append(np.double(content[1]))
        ref_plane.append(np.double(content[2]))
        ref_plane_list.append(ref_plane)
        
        line_id += 2
        n_coeff = int(lines[line_id])
        n_coeff_list.append(n_coeff)
        
        line_id += 2
        u_coeff = np.zeros((n_coeff,4))
        for i in range(n_coeff):
            u_coeff[i,:] = np.array(lines[line_id+i].split(',')).astype(np.double)
        u_coeff_list.append(u_coeff)
        
        v_coeff = np.zeros((n_coeff,4))
        for i in range(n_coeff):
            v_coeff[i,:] = np.array(lines[line_id+n_coeff+1+i].split(',')).astype(np.double)
        v_coeff_list.append(v_coeff)
#%%
digit = "{:.8e}"
for cam_id in range(ncam):
    file = '../../solutions/test_Camera/'+'cam'+str(cam_id+1)+'_poly'+'.txt'
    with open(file, 'w') as f:
        f.write('# Camera Model: (PINHOLE/POLYNOMIAL)\n' + str('POLYNOMIAL') + '\n')
        f.write('# Camera Calibration Error: \n' + str(None) + '\n')
        
        f.write('# Image Size: (n_row,n_col)\n')
        f.write(str(imgSize_list[cam_id][0]) + ',' + str(imgSize_list[cam_id][1]) + '\n')
        
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
