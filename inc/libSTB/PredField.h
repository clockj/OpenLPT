#ifndef PREDFIELD_H
#define PREDFIELD_H

#include <time.h>
#include <vector>
#include <algorithm>
#include <omp.h>

#include "Matrix.h"
#include "myMATH.h"
#include "STBCommons.h"
#include "ObjectInfo.h"


struct PFParam
{
    // X,Y,Z limits of view area
    AxisLimit limit;
    // # of grids in X,Y,Z >=2
    int nx, ny, nz;
    
    // radius of search sphere
    //  usually set: r <= 0.5*min(dx,dy,dz)
    //  it is better to make sure displacement < 0.5 * r
    double r;
    
    // number of bins for displacement statistics
    int nBin_x = 10; // >=2, if <2, set to 2
    int nBin_y = 10; // >=2, if <2, set to 2
    int nBin_z = 10; // >=2, if <2, set to 2
};



class PredField
{
public:
    PFParam _param;
    Matrix<double> _disp_field; // _nGrid_tot*3

    template<class T3D>
    PredField (PFParam const& param, std::vector<T3D> const& obj3d_list_prev, std::vector<T3D> const& obj3d_list_curr);
    PredField (PFParam const& param, std::vector<Pt3D> const& pt3d_list_prev, std::vector<Pt3D> const& pt3d_list_curr);

    // Directly set displacement field
    PredField (PFParam const& param, Matrix<double> const& disp_field);

    ~PredField() {};

    // interpolate displacement at (x,y,z)
    void getDisp (Pt3D& disp, Pt3D const& pt3d);

private:
    int _nGrid_tot; // total number of grid points
    double _dx, _dy, _dz; // grid spacing in X,Y,Z
    std::vector<double> _grid_x;
    std::vector<double> _grid_y;
    std::vector<double> _grid_z;

    // Calculate displacement (dx,dy,dz) statistics at each grid point
    // dx,dy,dz in [-2r,2r]
    // bin_center = -2r + 4r*bin_id/(disp_bins-1), bin_id in [0,disp_bins-1]
    //  e.g. n_bins = 3, bin_center = [-2r,0,2r], bin_id = [0,1,2]
    //  _size = n_bins: # of bins >=2
    //  _m = (_size-1)/(4*r): 1/bin_width
    //  _c = (_size-1)/2
    //  bin_id = round(_m*disp+_c)
    int _nBin_tot; // # of bins >=2
    std::vector<double> _m_xyz; // 1/bin_width
    std::vector<double> _c_xyz; // (_size-1)/2 

    // set grid points
    void setGrid();

    // create displacement field
    template<class T3D>
    void calDispField (std::vector<T3D> const& obj3d_list_prev, std::vector<T3D> const& obj3d_list_curr);
    template<class T3D>
    void calDispGrid (std::vector<int> const& xyz_id, std::vector<T3D> const& obj3d_list_prev, std::vector<T3D> const& obj3d_list_curr);

    void calDispField (std::vector<Pt3D> const& pt3d_list_prev, std::vector<Pt3D> const& pt3d_list_curr);
    void calDispGrid (std::vector<int> const& xyz_id, std::vector<Pt3D> const& pt3d_list_prev, std::vector<Pt3D> const& pt3d_list_curr);

    // find pt id within a sphere of radius r
    template<class T3D>
    void findNeighbor (std::vector<int>& pt_list_id, std::vector<int> const& xyz_id, std::vector<T3D> const& obj3d_list);
    void findNeighbor (std::vector<int>& pt_list_id, std::vector<int> const& xyz_id, std::vector<Pt3D> const& pt3d_list);

    // map from x_id,y_id,z_id to grid_id 
    // i in 0 ~ n_grid: i = index_x*_ny*_nz + index_y*_nz + index_z
    // for index_x in 0:_nx-1 
    //     for index_y in 0:_ny-1 
    //         for index_z in 0:_nz-1  
    //             i ++
    int mapGridID (int x_id, int y_id, int z_id);

    // map from xBin_id, yBin_id, zBin_id to bin_id
    // bin_id in 0 ~ nBin_tot: 
    // bin_id = xBin_id*nBin_y*nBin_z + yBin_id*nBin_z + zBin_id
    // for xBin_id in 0:nBin_x-1 
    //     for yBin_id in 0:nBin_y-1 
    //         for zBin_id in 0:nBin_z-1  
    //             bin_id ++    
    int mapBinID (int xBin_id, int yBin_id, int zBin_id);

    // update the "pdf" of each bin by using gaussian kernel
    // dist_pdf: used to find the most probable displacement, no need to be normalized
    void updateDispPDF (std::vector<double>& disp_pdf, std::vector<std::vector<double>> const& disp_list);

    // find the peak location of displacement pdf
    void findPDFPeakLoc (std::vector<double>& dist_opt, std::vector<double> const& disp_pdf);

    // find the peak bin id of displacement pdf
    void findPDFPeakID (std::vector<int>& peak_id, std::vector<double> const& disp_pdf);

    // fit the peak location based on guassian distribution
    double fitPeakBinLocGauss (double y1, double v1, double y2, double v2, double y3, double v3);

};

#include "PredField.hpp"

#endif