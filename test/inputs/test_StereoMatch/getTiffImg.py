#%%
import numpy as np
import cv2
from scipy.io import loadmat
import matplotlib.pyplot as plt

# %%
    
def getTiffImg(pt3d_list, rotVec, transVec, camMat, distCoeff, imgNRowNCol):
    npts = pt3d_list.shape[0]
    
    pt2d_list = cv2.projectPoints(pt3d_list.reshape(npts,1,3), rotVec, transVec, camMat, distCoeff)[0].reshape(npts, 2)
    
    # tracer radius
    tr_radius = 2

    # gaussian intensity distribution parameter
    alpha = 0
    a = 125
    b = 1.5
    c = 1.5
    min_intensity = 0
    max_intensity = 255
    dtype = np.uint8
    
    # calculate intensity
    img = np.zeros((imgNRowNCol[0], imgNRowNCol[1]), dtype=np.double)
    for i in range(npts):
        x = pt2d_list[i,0]
        y = pt2d_list[i,1]
        
        xmin = int(np.floor(max(x - tr_radius, 0)))
        xmax = int(np.floor(min(x + tr_radius + 1, imgNRowNCol[1])))
        ymin = int(np.floor(max(y - tr_radius, 0)))
        ymax = int(np.floor(min(y + tr_radius + 1, imgNRowNCol[0])))
        
        for j in range(ymin, ymax):
            for k in range(xmin, xmax):
                kk = (k-x) * np.cos(alpha) + (j-y) * np.sin(alpha)
                jj = -(k-x) * np.sin(alpha) + (j-y) * np.cos(alpha)
                
                img[j,k] = max(img[j,k], a * np.exp(-b*jj*jj-c*kk*kk))
                img[j,k] = np.clip(img[j,k], min_intensity, max_intensity)
    
    return img.astype(dtype)