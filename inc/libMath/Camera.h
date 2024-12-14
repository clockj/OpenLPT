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
#include <tuple>

#include "myMATH.h"
#include "Matrix.h"
#include "STBCommons.h"

struct PinholeParam
{
    int n_row; // number of rows of the image
    int n_col; // number of columns of the image
    Matrix<double> cam_mtx; // camera matrix (intrinsic parameters)
    bool is_distorted; // whether the image is distorted
    int n_dist_coeff; // number of distortion coefficients (4,5,8,12)
    std::vector<double> dist_coeff; // distortion coefficients
    Matrix<double> r_mtx;   // R matrix, rotation matrix (world coordinate -> camera coordinate)       
    Pt3D t_vec;   // T vector, translation vector (world coordinate -> camera coordinate)
    Matrix<double> r_mtx_inv; // inverse(R)
    Pt3D t_vec_inv; // - inv(R) @ T: center of camera in world coordinate   
};

enum RefPlane
{
    REF_X = 1,
    REF_Y,
    REF_Z
};
struct PolyParam
{
    int n_row; // number of rows of the image
    int n_col; // number of columns of the image
    RefPlane ref_plane; // reference plane
    std::vector<double> plane = {0,0}; // plane locations
    int n_coeff; // number of coefficients
    Matrix<double> u_coeffs; // coeff,x_power, y_power, z_power
    Matrix<double> du_coeffs; // du/dvar1 coeff,x_power, y_power, z_power
                              // du/dvar2 coeff,x_power, y_power, z_power
    Matrix<double> v_coeffs; // coeff,x_power, y_power, z_power 
    Matrix<double> dv_coeffs; // dv/dvar1 coeff,x_power, y_power, z_power
                              // dv/dvar2 coeff,x_power, y_power, z_power
};

struct PinPlateParam: public PinholeParam
{ 
    Plane3D plane; // reference refractive plane (farthest plane to camera), normal vector is pointing away from the camera
    Pt3D pt3d_closest; // closest point to the camera
    std::vector<double> refract_array; // refractive index array (from farthest to nearest)
    std::vector<double> w_array; // width of the refractive plate (from farthest to nearest)
    int n_plate; // number of refractive plates
    double proj_tol; // projection tolerance squatre
    int proj_nmax; // maximum number of iterations for projection
    double lr; // learning rate
    double refract_ratio_max; // max(refract_array[0]/[i])
};

enum CameraType
{
    PINHOLE = 0,
    POLYNOMIAL,
    PINPLATE
};

class Camera
{
public:
    CameraType   _type;
    PinholeParam _pinhole_param;
    PolyParam    _poly_param; 
    PinPlateParam _pinplate_param;

    Camera ();
    Camera (const Camera& c); // Camera deep copy
    Camera (std::istream& is);
    Camera (std::string file_name);
    ~Camera () {};

    // Load parameters from an input string 
    // input: is (file Input stream) 
    void loadParameters (std::istream& is);
    void loadParameters (std::string file_name);

    void updatePolyDuDv ();

    // Transfer from rotation matrix to rotation vector
    Pt3D rmtxTorvec (Matrix<double> const& r_mtx);

    // Save parameters to an output string
    void saveParameters (std::string file_name);


    //                //
    // Get image size //
    //                //
    int getNRow () const; // get number of rows of the image
    int getNCol () const; // get number of columns of the image


    //            //
    // Projection //
    //            //
    Pt2D project (Pt3D const& pt_world, bool is_print_detail = false) const;

    // Project world coordinate [mm] to image coordinate [mm]: 
    //  (xw,yw,zw) -> (x,y,z) -> (Xu,Yu,0)
    //  (x,y,z) = R @ (xw,yw,zw) + T
    //  (Xu,Yu,0) = (fx*x/z + cx, fy*y/z + cy)
    // input: pt_world: point location in world coordinate 
    // output
    Pt2D worldToUndistImg (Pt3D const& pt_world, PinholeParam const& param) const;

    std::tuple<bool, Pt3D, double> refractPlate (Pt3D const& pt_world) const;

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
    Pt2D distort (Pt2D const& pt_img_undist, PinholeParam const& param) const;

    // Polynomial model 
    // Project world coordinate [mm] to distorted image [px]:
    //  (xw,yw,zw) -> (xd,yd)
    Pt2D polyProject (Pt3D const& pt_world) const;


    //               //
    // Line of sight //
    //               //
    Line3D lineOfSight (Pt2D const& pt_img_undist) const;

    // Project Image pixel to camera coordinate [mm]: 
    //  (Xf,Yf,0) -> (Xd,Yd,0) -> (Xu,Yu,0)
    // input: pt_pix: point location in pixel unit on image 
    //        (no use of z coordinate of pt_pix)
    // output: Matrix (camera coordinate) (Xu,Yu,0)
    Pt2D undistort (Pt2D const& pt_img_dist, PinholeParam const& param) const;

    // Project image coordinate [mm] to world coordinate [mm]: 
    //  (Xu,Yu,0) -> (x,y,z) -> (xw,yw,zw)
    // input: pt_img_undist: undistorted image coordinate 
    // output: line of sight
    Line3D pinholeLine (Pt2D const& pt_img_undist) const;

    Line3D pinplateLine (Pt2D const& pt_img_undist) const;

    // Polynomial Model
    Pt3D polyImgToWorld (Pt2D const& pt_img_dist, double plane_world) const;

    Line3D polyLineOfSight (Pt2D const& pt_img_dist) const;

};

struct CamList
{
    std::vector<Camera> cam_list; // cam_id: 0,1,2,3,..
    std::vector<int> intensity_max; // default: 8 digit (255), save size as cam_list
    std::vector<int> useid_list; // cam_id to be used
};

#endif