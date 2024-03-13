#ifndef OTF_H
#define OTF_H

#include <vector>
#include "Matrix.h"
#include "STBCommons.h"
#include "myMATH.h"

struct OTFParam
{
    double dx, dy, dz;
    int n_cam, nx, ny, nz, n_grid; // n_grid = nx*ny*nz

    // OTF parameters: calculate intensity of a tracer onto camera
    //  2D projection point: (x0, y0)
    //  xx =  (x-x0) * cos(alpha) + (y-y0) * sin(alpha)
    //  yy = -(x-x0) * sin(alpha) + (y-y0) * cos(alpha)
    //  img[y,x] = a * exp(-b * xx^2 - c * yy^2) 
    // default: a = 125, b = 1.5, c = 1.5, alpha = 0   
    Matrix<double> a, b, c, alpha; 
    
    AxisLimit boundary;
    std::vector<double> grid_x, grid_y, grid_z;
};

class OTF 
{
private:
    //  i in 0 ~ n_grid: i = x_id * (ny*nz) + y_id * nz + z_id
    //  for x_id = 0:_nx-1 
    //      for y_id = 0:_ny-1 
    //          for z_id = 0:_nz-1  
    //              i ++
    int mapGridID (int x_id, int y_id, int z_id) const
    {
        int id = x_id * (_param.ny*_param.nz) + y_id * _param.nz + z_id;
        return id;
    }
    void setGrid ();

public:
    OTFParam _param;

    OTF () {};
    OTF(const OTF& otf) : _param(otf._param) {};

    // nx,ny,nz >= 2
    // set as default a,b,c,alpha values
    OTF(int n_cam, int nx, int ny, int nz, AxisLimit const& boundary); 

    // nx,ny,nz >= 2
    // use file to load a,b,c,alpha
    OTF(std::string otf_file); 

    ~OTF() {};

    void loadParam (int n_cam, int nx, int ny, int nz, AxisLimit const& boundary);
    void loadParam (std::string otf_file);

    // Output: (a,b,c,alpha)
    std::vector<double> getOTFParam(int cam_id, Pt3D const& pt_world) const;

    // TODO: add self-calibration code.
};

#endif
