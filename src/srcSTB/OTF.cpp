#include "OTF.h"

OTF::OTF(int n_cam, int nx, int ny, int nz, AxisLimit const& boundary)
{
    loadParam(n_cam, nx, ny, nz, boundary);
};

OTF::OTF(std::string otf_file)
{
    loadParam(otf_file);
};

void OTF::setGrid ()
{
    _param.grid_x = myMATH::linspace(_param.boundary.x_min, _param.boundary.x_max, _param.nx);
    _param.grid_y = myMATH::linspace(_param.boundary.y_min, _param.boundary.y_max, _param.ny);
    _param.grid_z = myMATH::linspace(_param.boundary.z_min, _param.boundary.z_max, _param.nz);

    _param.dx = (_param.boundary.x_max - _param.boundary.x_min) / (_param.nx - 1);
    _param.dy = (_param.boundary.y_max - _param.boundary.y_min) / (_param.ny - 1);
    _param.dz = (_param.boundary.z_max - _param.boundary.z_min) / (_param.nz - 1);
}

void OTF::loadParam (int n_cam, int nx, int ny, int nz, AxisLimit const& boundary)
{
    if (nx < 2 || ny < 2 || nz < 2)
    {
        std::cerr << "OTF::loadParam error at line" << __LINE__ << ":\n"
                  << "nx, ny, nz should be larger than 1." << std::endl;
        std::cerr << "nx = " << nx << ", ny = " << ny << ", nz = " << nz << std::endl;
        throw error_range;
    }

    _param.n_cam = n_cam;
    _param.nx = nx;
    _param.ny = ny;
    _param.nz = nz;
    _param.n_grid = nx * ny * nz;

    _param.boundary = boundary;

    _param.a = Matrix<double> (n_cam, _param.n_grid, 125);
    _param.b = Matrix<double> (n_cam, _param.n_grid, 1.5);
    _param.c = Matrix<double> (n_cam, _param.n_grid, 1.5);
    _param.alpha = Matrix<double> (n_cam, _param.n_grid, 0);

    setGrid();
}

void OTF::loadParam (std::string otf_file)
{
    std::ifstream infile(otf_file);

    if (!infile.is_open())
    {
        std::cerr << "OTF::loadParam error at line" << __LINE__ << ":\n"
                  << "Cannot open file " << otf_file << std::endl;
        throw error_io;
    }

    std::string line;
    std::stringstream file_content;
    while (std::getline(infile, line))
    {
        size_t comment_pos = line.find("#");
        if (comment_pos > 0 && comment_pos < std::string::npos)
        {
            line.erase(comment_pos);
        }
        else if (comment_pos == 0)
        {
            continue;
        }

        file_content << line << '\t';
    }
    infile.close();
    
    file_content >> _param.n_cam;
    file_content.ignore();
    file_content >> _param.nx;
    file_content.ignore();
    file_content >> _param.ny;
    file_content.ignore();
    file_content >> _param.nz;
    file_content.ignore();
    file_content >> _param.n_grid;

    file_content >> _param.boundary.x_min;
    file_content.ignore();
    file_content >> _param.boundary.x_max;
    file_content.ignore();
    file_content >> _param.boundary.y_min;
    file_content.ignore();
    file_content >> _param.boundary.y_max;
    file_content.ignore();
    file_content >> _param.boundary.z_min;
    file_content.ignore();
    file_content >> _param.boundary.z_max;

    _param.a = Matrix<double> (_param.n_cam, _param.n_grid, file_content);
    _param.b = Matrix<double> (_param.n_cam, _param.n_grid, file_content);
    _param.c = Matrix<double> (_param.n_cam, _param.n_grid, file_content);
    _param.alpha = Matrix<double> (_param.n_cam, _param.n_grid, file_content);

    setGrid();
}

std::vector<double> OTF::getOTFParam(int cam_id, Pt3D const& pt3d) const
{ 
    double pt3d_x = std::min(
        std::max(pt3d[0], _param.boundary.x_min),
        _param.boundary.x_max
    );
    double pt3d_y = std::min(
        std::max(pt3d[1], _param.boundary.y_min),
        _param.boundary.y_max
    );
    double pt3d_z = std::min(
        std::max(pt3d[1], _param.boundary.z_min),
        _param.boundary.z_max
    );
    
    // find out the limits of interpolation cube
    int x_id = std::max(
        0, 
        (int) std::floor((pt3d_x - _param.boundary.x_min) / _param.dx)
    );
    int y_id = std::max(
        0, 
        (int) std::floor((pt3d_y - _param.boundary.y_min) / _param.dy)
    );
    int z_id = std::max(
        0, 
        (int) std::floor((pt3d_z - _param.boundary.z_min) / _param.dz)
    );
    x_id = std::min(x_id, _param.nx-2);
    y_id = std::min(y_id, _param.ny-2);
    z_id = std::min(z_id, _param.nz-2);

    AxisLimit grid_limit;

    grid_limit.x_min = _param.grid_x[x_id];
    grid_limit.x_max = _param.grid_x[x_id + 1];
    grid_limit.y_min = _param.grid_y[y_id];
    grid_limit.y_max = _param.grid_y[y_id + 1];
    grid_limit.z_min = _param.grid_z[z_id];
    grid_limit.z_max = _param.grid_z[z_id + 1];

    int i_000 = mapGridID(x_id  , y_id  , z_id  );
    int i_100 = mapGridID(x_id+1, y_id  , z_id  );
    int i_101 = mapGridID(x_id+1, y_id  , z_id+1);
    int i_001 = mapGridID(x_id  , y_id  , z_id+1);
    int i_010 = mapGridID(x_id  , y_id+1, z_id  );
    int i_110 = mapGridID(x_id+1, y_id+1, z_id  );
    int i_111 = mapGridID(x_id+1, y_id+1, z_id+1);
    int i_011 = mapGridID(x_id  , y_id+1, z_id+1);

    std::vector<double> a_value = { 
        _param.a(cam_id,i_000), 
        _param.a(cam_id,i_100),
        _param.a(cam_id,i_101),
        _param.a(cam_id,i_001),
        _param.a(cam_id,i_010),
        _param.a(cam_id,i_110),
        _param.a(cam_id,i_111),
        _param.a(cam_id,i_011)
    };
    std::vector<double> b_value = { 
        _param.b(cam_id,i_000), 
        _param.b(cam_id,i_100),
        _param.b(cam_id,i_101),
        _param.b(cam_id,i_001),
        _param.b(cam_id,i_010),
        _param.b(cam_id,i_110),
        _param.b(cam_id,i_111),
        _param.b(cam_id,i_011)
    };
    std::vector<double> c_value = { 
        _param.c(cam_id,i_000), 
        _param.c(cam_id,i_100),
        _param.c(cam_id,i_101),
        _param.c(cam_id,i_001),
        _param.c(cam_id,i_010),
        _param.c(cam_id,i_110),
        _param.c(cam_id,i_111),
        _param.c(cam_id,i_011)
    };
    std::vector<double> alpha_value = { 
        _param.alpha(cam_id,i_000), 
        _param.alpha(cam_id,i_100),
        _param.alpha(cam_id,i_101),
        _param.alpha(cam_id,i_001),
        _param.alpha(cam_id,i_010),
        _param.alpha(cam_id,i_110),
        _param.alpha(cam_id,i_111),
        _param.alpha(cam_id,i_011)
    };

    std::vector<double> res(4,0); // a,b,c,alpha
    std::vector<double> pt_vec = {pt3d_x, pt3d_y, pt3d_z};
    res[0] = myMATH::triLinearInterp(grid_limit, a_value, pt_vec);
    res[1] = myMATH::triLinearInterp(grid_limit, b_value, pt_vec);
    res[2] = myMATH::triLinearInterp(grid_limit, c_value, pt_vec);
    res[3] = myMATH::triLinearInterp(grid_limit, alpha_value, pt_vec);

    return res;
}