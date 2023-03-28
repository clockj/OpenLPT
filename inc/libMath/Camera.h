//
//  Camera.h
//  
//  header file of class Camera 
//  
//  Created  by Shijie Zhong on 07/02/2022
//

#ifndef CAMERA_H
#define CAMERA_H

#include <fstream>
#include <string.h>
#include <iostream>
#include <sstream>

#include "Matrix.h"
#include "myMATH.h"
#include "STBCommons.h"


class Camera
{
protected:
    // Model parameters; names should be more or less the same
    //  as from calibTsai.m! 
    double  _n_off_h = 0; // pixel center coordinate (Yf) 
    double  _n_off_w = 0; // pixel center coordinate (Xf)
    int     _n_pix_w = 0; // number of pixel along X direction (col)
    int     _n_pix_h = 0; // number of pixel along Y direction (row)
    double  _w_pix = 0;   // distance between 2 pixel centers along X direction [mm]
    double  _h_pix = 0;	  // distance between 2 pixel centers along Y direction [mm] 
    double  _f_eff = 0;   // effective focus length [mm]
    double  _kr = 0;      // distortion coefficient along radial direction
    double  _kx = 0;      // distortion coefficient along X direction
    Matrix<double> _r_mtx; // R matrix, rotation matrix (world coordinate -> camera coordinate)       
    Matrix<double> _t_vec; // T vector, translation vector (world coordinate -> camera coordinate)
    Matrix<double> _r_mtx_inv; // inverse(R)
    Matrix<double> _t_vec_inv; // inverse(T)

public:
    Camera ();
    Camera (std::istream& is);
    Camera (std::string file_name);
    // Camera deep copy
    Camera (const Camera& c); 
    ~Camera () {};

    // Load parameters from an input string 
    // input: is (file Input stream) 
    void LoadParameters (std::istream& is);
    void LoadParametersFromFile (std::string file_name);

    //            //
    // Projection //
    //            //

    // Project Image pixel to camera coordinate [mm]: 
    //  (Xf,Yf,0) -> (Xd,Yd,0) -> (Xu,Yu,0)
    // input: pt_pix: point location in pixel unit on image 
    //        (no use of z coordinate of pt_pix)
    // output: Matrix (camera coordinate) (Xu,Yu,0)
    Matrix<double> ImgPixelToMM (Matrix<double> pt_img_pix);

    // Project Image in camera coordinate [mm] to pixel: 
    //  (Xu,Yu,0) -> (Xd,Yd,0) -> (Xf,Yf,0)
    // input: pt_cam: point location in camera coordinate on image 
    //        (no use of z coordinate of pt_cam)
    // output: Matrix (camera coordinate) (Xf,Yf,0)
    // origin is on the top left corner
    // ----> x direction (rightwards)
    // | 
    // |
    // \/ y direction (downwards)
    Matrix<double> ImgMMToPixel (Matrix<double> pt_img_mm);

    // Project image coordinate [mm] to world coordinate [mm]: 
    //  (Xu,Yu,0) -> (x,y,z) -> (xw,yw,zw)
    // input: pt_cam: point location in world coordinate 
    //        (no use of z coordinate of pt_cam)
    // output: Matrix (camera coordinate) (xw,yw,zw)
    Matrix<double> ImgMMToWorld (Matrix<double> pt_img_mm);

    // Project world coordinate [mm] to image coordinate [mm]: 
    //  (xw,yw,zw) -> (x,y,z) -> (Xu,Yu,0)
    // input: pt_world: point location in world coordinate 
    //        (no use of z coordinate of pt_cam)
    // output: Matrix (camera coordinate) (xw,yw,zw)
    Matrix<double> WorldToImgMM (Matrix<double> pt_world);



    //
    // Get the parameters of the camera images 
    //

    double GetNoffh () {return _n_off_h;};
    double GetNoffw () {return _n_off_w;};
    int    GetNpixw () {return _n_pix_w;};
    int    GetNpixh () {return _n_pix_h;};
    double GetWpix  () {return _w_pix;};
    double GetHpix  () {return _h_pix;};
    double GetKr    () {return _kr;};
    double GetKx    () {return _kx;};

    // Return camera projective center (in world coordinates)
    // output: Matrix (xw;yw;zw)
    Matrix<double> GetCenterInWorld () {return _t_vec_inv;};
};

# endif