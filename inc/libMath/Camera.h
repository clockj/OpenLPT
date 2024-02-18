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
#include <string>
#include <iostream>
#include <sstream>
#include <limits>

#include "myMATH.h"
#include "Matrix.h"
#include "STBCommons.h"

struct PinholeParam
{
    Matrix<double> cam_mtx; // camera matrix (intrinsic parameters)
    bool is_distorted; // whether the image is distorted
    int n_dist_coeff; // number of distortion coefficients (4,5,8,12)
    std::vector<double> dist_coeff; // distortion coefficients
    Matrix<double> r_mtx;   // R matrix, rotation matrix (world coordinate -> camera coordinate)       
    Pt3D t_vec;   // T vector, translation vector (world coordinate -> camera coordinate)
    Matrix<double> r_mtx_inv; // inverse(R)
    Pt3D t_vec_inv; // - inv(R) @ T: center of camera in world coordinate   
};

struct PolyParam
{
    std::vector<double> coeff;
    Matrix<double> order; // order of polynomial in x, y, z
};

enum CameraType
{
    PINHOLE = 0,
    POLYNOMIAL
};

class Camera
{
public:
    CameraType   _type;
    PinholeParam _pinhole_param;
    PolyParam    _poly_param; 

    Camera ();
    Camera (const Camera& c); // Camera deep copy
    Camera (std::istream& is);
    Camera (std::string file_name);
    ~Camera () {};

    // Load parameters from an input string 
    // input: is (file Input stream) 
    void loadParameters (std::istream& is);
    void loadParameters (std::string file_name);

    // Transfer from rotation matrix to rotation vector
    Pt3D rmtxTorvec (Matrix<double> const& r_mtx);

    // Save parameters to an output string
    void saveParameters (std::string file_name);

    //            //
    // Projection //
    //            //
    Pt2D project (Pt3D const& pt_world);

    // Project world coordinate [mm] to image coordinate [mm]: 
    //  (xw,yw,zw) -> (x,y,z) -> (Xu,Yu,0)
    //  (x,y,z) = R @ (xw,yw,zw) + T
    //  (Xu,Yu,0) = (fx*x/z + cx, fy*y/z + cy)
    // input: pt_world: point location in world coordinate 
    // output
    Pt2D worldToUndistImg (Pt3D const& pt_world);

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
    Pt2D distort (Pt2D const& pt_img_undist);


    //               //
    // Line of sight //
    //               //
    Line3D lineOfSight (Pt2D const& pt_img_undist);

    // Project Image pixel to camera coordinate [mm]: 
    //  (Xf,Yf,0) -> (Xd,Yd,0) -> (Xu,Yu,0)
    // input: pt_pix: point location in pixel unit on image 
    //        (no use of z coordinate of pt_pix)
    // output: Matrix (camera coordinate) (Xu,Yu,0)
    Pt2D undistort (Pt2D const& pt_img_dist);

    // Project image coordinate [mm] to world coordinate [mm]: 
    //  (Xu,Yu,0) -> (x,y,z) -> (xw,yw,zw)
    // input: pt_img_undist: undistorted image coordinate 
    // output: line of sight
    Line3D pinholeLine (Pt2D const& pt_img_undist);



};

# endif