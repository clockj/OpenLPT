#ifndef OTF_H
#define OTF_H

#include <vector>
#include "Matrix.h"
#include "STBCommons.h"
#include "myMATH.h"

class OTF 
{
protected:
    // _a, _b, _c, _alpha: n_cam*n_grid, n_grid = _nx*_ny*_nz
    //  i in 0 ~ n_grid: i = index_x*_ny*_nz + index_y*_nz + index_z
    //  for index_x in 0:_nx 
    //      for index_y in 0:_ny 
    //          for index_z in 0:_nz  
    //              i ++
    int _n_cam, _nx, _ny, _nz;
    Matrix<double> _a, _b, _c, _alpha;
    double _dx, _dy, _dz;
    std::vector<double> _grid_x, _grid_y, _grid_z;
    AxisLimit _boundary;

    void InitParam (double a_val, double b_val, double c_val, double alpha_val);

public:
    // TODO: add self-calibration code.
    // TODO: run more tests!

    // OTF() {};
    ~OTF() {};

    OTF(int n_cam, int n, AxisLimit& boundary) 
        : _n_cam(n_cam), _nx(n), _ny(n), _nz(n) 
    {
        _boundary = boundary;
        InitParam(125, 1.5, 1.5, 0);
    };
    
    OTF(int n_cam, int n_x, int n_y, int n_z, AxisLimit& boundary) 
        : _n_cam(n_cam), _nx(n_x), _ny(n_y), _nz(n_z) 
    {
        _boundary = boundary;
        InitParam(125, 1.5, 1.5, 0);
    };

    OTF(const OTF& otf)
        : _a(otf._a), _b(otf._b), _c(otf._c), _alpha(otf._alpha), 
          _n_cam(otf._n_cam), _nx(otf._nx), _ny(otf._ny), _nz(otf._nz),
          _dx(otf._dx), _dy(otf._dy), _dz(otf._dz), 
          _grid_x(otf._grid_x), _grid_y(otf._grid_y), _grid_z(otf._grid_z)
    {
        _boundary = otf._boundary;
    };

    int MapGridIndexToID (int index_x, int index_y, int index_z)
    {
        return index_x*_ny*_nz + index_y*_nz + index_z;
    }

    // Output: (a,b,c,alpha)
    std::vector<double> GetOTFParam(int cam_index, const Matrix<double>& pt_world);

};

#endif
