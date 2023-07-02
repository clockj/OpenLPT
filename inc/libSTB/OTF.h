#ifndef OTF_H
#define OTF_H

#include <vector>
#include "Matrix.h"
#include "STBCommons.h"
#include "myMATH.h"

class OTF 
{
private:
    // _a, _b, _c, _alpha: n_cam*n_grid, n_grid = _nx*_ny*_nz
    //  i in 0 ~ n_grid: i = index_x*_ny*_nz + index_y*_nz + index_z
    //  for index_x in 0:_nx-1 
    //      for index_y in 0:_ny-1 
    //          for index_z in 0:_nz-1  
    //              i ++
    int _n_cam, _nx, _ny, _nz, _n_grid;
    Matrix<double> _a, _b, _c, _alpha;
    double _dx, _dy, _dz;
    AxisLimit _boundary;
    std::vector<double> _grid_x, _grid_y, _grid_z;

    void LoadParam (int n_cam, int n_x, int n_y, int n_z, AxisLimit& boundary);
    void LoadParam (int n_cam, int n, AxisLimit& boundary);
    void SetGrid ();

public:
    // TODO: add self-calibration code.
    OTF () {};
    OTF(int n_cam, int n_x, int n_y, int n_z, AxisLimit& boundary) 
        : _n_cam(n_cam), _nx(n_x), _ny(n_y), _nz(n_z),
          _n_grid(_nx*_ny*_nz),
          _a(_n_cam, _n_grid, 125),
          _b(_n_cam, _n_grid, 1.5),
          _c(_n_cam, _n_grid, 1.5),
          _alpha(_n_cam, _n_grid, 0),
          _boundary(boundary)
    {
        SetGrid();
    };
    OTF(int n_cam, int n, AxisLimit& boundary) 
        : OTF(n_cam, n, n, n, boundary)
    {};

    OTF(int n_cam, int n_x, int n_y, int n_z, AxisLimit& boundary, std::string file_head);
    OTF(int n_cam, int n, AxisLimit& boundary, std::string file_head)
        : OTF(n_cam, n, n, n, boundary, file_head)
    {};

    OTF(const OTF& otf)
        : _n_cam(otf._n_cam), _nx(otf._nx), _ny(otf._ny), _nz(otf._nz),
          _n_grid(otf._n_grid),_a(otf._a), _b(otf._b), _c(otf._c), _alpha(otf._alpha), 
          _dx(otf._dx), _dy(otf._dy), _dz(otf._dz), _boundary(otf._boundary),
          _grid_x(otf._grid_x), _grid_y(otf._grid_y), _grid_z(otf._grid_z)
    {};

    ~OTF() {};

    int MapGridIndexToID (int index_x, int index_y, int index_z)
    {
        return index_x*_ny*_nz + index_y*_nz + index_z;
    }

    // Output: (a,b,c,alpha)
    std::vector<double> GetOTFParam(int cam_index, const Matrix<double>& pt_world);
    Matrix<double> GetA() 
    {
        return _a;
    }

};

#endif
