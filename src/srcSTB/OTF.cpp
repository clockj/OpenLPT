#include "OTF.h"


OTF::OTF(int n_cam, int n_x, int n_y, int n_z, AxisLimit& boundary, std::string file_head)
        : _n_cam(n_cam), _nx(n_x), _ny(n_y), _nz(n_z),
          _n_grid(_nx*_ny*_nz), 
          _a(file_head+"_a.csv"),_b(file_head+"_b.csv"),_c(file_head+"_c.csv"),_alpha(file_head+"_alpha.csv"),
          _boundary(boundary)
{
    bool judge = _a.GetDimX()!=_n_cam || _a.GetDimY()!=_n_grid || 
                 _b.GetDimX()!=_n_cam || _b.GetDimY()!=_n_grid ||
                 _c.GetDimX()!=_n_cam || _c.GetDimY()!=_n_grid ||
                 _alpha.GetDimX()!=_n_cam || _alpha.GetDimY()!=_n_grid;
    if (judge)
    {
        std::cerr << "OTF::OTF error: a/b/c/alpha dim is inconsistent!" << std::endl;
        std::cerr << "_n_cam,_n_grid=" << _n_cam << "," << _n_grid << std::endl;
        std::cerr << _a.GetDimX() << "," << _a.GetDimY() << ","
                    << _b.GetDimX() << "," << _b.GetDimY() << ","
                    << _c.GetDimX() << "," << _c.GetDimY() << ","
                    << _alpha.GetDimX() << "," << _alpha.GetDimY()
                    << std::endl;
        throw error_size;
    }

    SetGrid();
};


void OTF::LoadParam (int n_cam, int n_x, int n_y, int n_z, AxisLimit& boundary)
{
    _n_cam = n_cam; 
    _nx = n_x;
    _ny = n_y;
    _nz = n_z;
    _n_grid = _nx*_ny*_nz;
    _a = Matrix<double>(_n_cam, _n_grid, 125);
    _b = Matrix<double>(_n_cam, _n_grid, 1.5);
    _c = Matrix<double>(_n_cam, _n_grid, 1.5);
    _alpha = Matrix<double>(_n_cam, _n_grid, 0);
    _boundary = boundary;

    SetGrid();
}


void OTF::LoadParam (int n_cam, int n, AxisLimit& boundary)
{
    LoadParam (n_cam, n, n, n, boundary);
}


void OTF::SetGrid ()
{
    _grid_x = myMATH::Linspace(_boundary._x_min, _boundary._x_max, _nx);
    _grid_y = myMATH::Linspace(_boundary._y_min, _boundary._y_max, _ny);
    _grid_z = myMATH::Linspace(_boundary._z_min, _boundary._z_max, _nz);

    _dx = (_boundary._x_max - _boundary._x_min) / (_nx - 1);
    _dy = (_boundary._y_max - _boundary._y_min) / (_ny - 1);
    _dz = (_boundary._z_max - _boundary._z_min) / (_nz - 1);
}


std::vector<double> OTF::GetOTFParam(int cam_index, const Matrix<double>& pt_world)
{
    // find out the limits of interpolation cube 
    double pt_x = std::min(
        std::max(pt_world(0,0), _boundary._x_min),
        _boundary._x_max
    );
    double pt_y = std::min(
        std::max(pt_world(1,0), _boundary._y_min),
        _boundary._y_max
    );
    double pt_z = std::min(
        std::max(pt_world(2,0), _boundary._z_min),
        _boundary._z_max
    );
    
    int index_x = std::max(
        0, 
        (int) std::floor((pt_x - _boundary._x_min) / _dx)
    );
    int index_y = std::max(
        0, 
        (int) std::floor((pt_y - _boundary._y_min) / _dy)
    );
    int index_z = std::max(
        0, 
        (int) std::floor((pt_z - _boundary._z_min) / _dz)
    );

    AxisLimit grid_limit;

    grid_limit._x_min = _grid_x[index_x];
    grid_limit._x_max = _grid_x[index_x + 1];

    grid_limit._y_min = _grid_y[index_y];
    grid_limit._y_max = _grid_y[index_y + 1];

    grid_limit._z_min = _grid_z[index_z];
    grid_limit._z_max = _grid_z[index_z + 1];

    int i_000 = MapGridIndexToID(index_x  , index_y  , index_z  );
    int i_100 = MapGridIndexToID(index_x+1, index_y  , index_z  );
    int i_101 = MapGridIndexToID(index_x+1, index_y  , index_z+1);
    int i_001 = MapGridIndexToID(index_x  , index_y  , index_z+1);
    int i_010 = MapGridIndexToID(index_x  , index_y+1, index_z  );
    int i_110 = MapGridIndexToID(index_x+1, index_y+1, index_z  );
    int i_111 = MapGridIndexToID(index_x+1, index_y+1, index_z+1);
    int i_011 = MapGridIndexToID(index_x  , index_y+1, index_z+1);

    std::vector<double> a_value = { 
        _a(cam_index,i_000), 
        _a(cam_index,i_100),
        _a(cam_index,i_101),
        _a(cam_index,i_001),
        _a(cam_index,i_010),
        _a(cam_index,i_110),
        _a(cam_index,i_111),
        _a(cam_index,i_011)
    };
    std::vector<double> b_value = { 
        _b(cam_index,i_000), 
        _b(cam_index,i_100),
        _b(cam_index,i_101),
        _b(cam_index,i_001),
        _b(cam_index,i_010),
        _b(cam_index,i_110),
        _b(cam_index,i_111),
        _b(cam_index,i_011)
    };
    std::vector<double> c_value = { 
        _c(cam_index,i_000), 
        _c(cam_index,i_100),
        _c(cam_index,i_101),
        _c(cam_index,i_001),
        _c(cam_index,i_010),
        _c(cam_index,i_110),
        _c(cam_index,i_111),
        _c(cam_index,i_011)
    };
    std::vector<double> alpha_value = { 
        _alpha(cam_index,i_000), 
        _alpha(cam_index,i_100),
        _alpha(cam_index,i_101),
        _alpha(cam_index,i_001),
        _alpha(cam_index,i_010),
        _alpha(cam_index,i_110),
        _alpha(cam_index,i_111),
        _alpha(cam_index,i_011)
    };

    std::vector<double> res(4,0);
    std::vector<double> pt_vec = {pt_x, pt_y, pt_z};
    res[0] = myMATH::TriLinearInterp(grid_limit, a_value, pt_vec);
    res[1] = myMATH::TriLinearInterp(grid_limit, b_value, pt_vec);
    res[2] = myMATH::TriLinearInterp(grid_limit, c_value, pt_vec);
    res[3] = myMATH::TriLinearInterp(grid_limit, alpha_value, pt_vec);

    return res;
}