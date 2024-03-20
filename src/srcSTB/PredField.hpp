#ifndef PREDFIELD_HPP
#define PREDFIELD_HPP

#include "PredField.h"

template<class T3D>
PredField::PredField (PFParam const& param, std::vector<T3D> const& obj3d_list_prev, std::vector<T3D> const& obj3d_list_curr) : _param(param)
{
    setGrid();

    calDispField(obj3d_list_prev, obj3d_list_curr);
}


// For Pt3D list
PredField::PredField (PFParam const& param, std::vector<Pt3D> const& pt3d_list_prev, std::vector<Pt3D> const& pt3d_list_curr) : _param(param)
{
    setGrid();

    calDispField(pt3d_list_prev, pt3d_list_curr);
}


// Directly set displacement field
PredField::PredField (PFParam const& param, Matrix<double> const& disp_field) : _param(param)
{
    setGrid();

    if (disp_field.getDimRow() != _nGrid_tot || disp_field.getDimCol() != 3)
    {
        std::cerr << "PredField::PredField error at line " << __LINE__ << ":\n"
                  << "The input displacement field has wrong dimension." << "\n"
                  << "Expected: (" << _nGrid_tot << ",3), "
                  << "Got: (" << disp_field.getDimRow() << "," << disp_field.getDimCol() << ")." << std::endl;
        throw error_size;
    }

    _disp_field = disp_field;
}


void PredField::getDisp (Pt3D& disp, Pt3D const& pt3d)
{
    // find out the limits of interpolation cube 
    double pt_x = std::min(
        std::max(pt3d[0], _param.limit.x_min),
        _param.limit.x_max
    );
    double pt_y = std::min(
        std::max(pt3d[1], _param.limit.y_min),
        _param.limit.y_max
    );
    double pt_z = std::min(
        std::max(pt3d[2], _param.limit.z_min),
        _param.limit.z_max
    );

    int x_id = std::max(
        0, 
        (int) std::floor((pt_x - _param.limit.x_min) / _dx)
    );
    int y_id = std::max(
        0, 
        (int) std::floor((pt_y - _param.limit.y_min) / _dy)
    );
    int z_id = std::max(
        0, 
        (int) std::floor((pt_z - _param.limit.z_min) / _dz)
    );

    x_id = std::min(x_id, _param.nx-2);
    y_id = std::min(y_id, _param.ny-2);
    z_id = std::min(z_id, _param.nz-2);

    AxisLimit grid_limit;

    grid_limit.x_min = _grid_x[x_id];
    grid_limit.x_max = _grid_x[x_id + 1];

    grid_limit.y_min = _grid_y[y_id];
    grid_limit.y_max = _grid_y[y_id + 1];

    grid_limit.z_min = _grid_z[z_id];
    grid_limit.z_max = _grid_z[z_id + 1];

    int i_000 = mapGridID(x_id  , y_id  , z_id  );
    int i_100 = mapGridID(x_id+1, y_id  , z_id  );
    int i_101 = mapGridID(x_id+1, y_id  , z_id+1);
    int i_001 = mapGridID(x_id  , y_id  , z_id+1);
    int i_010 = mapGridID(x_id  , y_id+1, z_id  );
    int i_110 = mapGridID(x_id+1, y_id+1, z_id  );
    int i_111 = mapGridID(x_id+1, y_id+1, z_id+1);
    int i_011 = mapGridID(x_id  , y_id+1, z_id+1);

    std::vector<double> pt_vec = {pt_x, pt_y, pt_z};

    for (int j = 0; j < 3; j ++)
    {
        std::vector<double> field = {
            _disp_field(i_000, j),
            _disp_field(i_100, j),
            _disp_field(i_101, j),
            _disp_field(i_001, j),
            _disp_field(i_010, j),
            _disp_field(i_110, j),
            _disp_field(i_111, j),
            _disp_field(i_011, j)
        };
        disp[j] = myMATH::triLinearInterp(grid_limit, field, pt_vec);
    }
}


void PredField::saveDispField (std::string const& file)
{
    _disp_field.write(file);
}


void PredField::setGrid()
{
    _nGrid_tot = _param.nx*_param.ny*_param.nz;
    _dx = (_param.limit.x_max - _param.limit.x_min)/(_param.nx-1);
    _dy = (_param.limit.y_max - _param.limit.y_min)/(_param.ny-1);
    _dz = (_param.limit.z_max - _param.limit.z_min)/(_param.nz-1);

    _grid_x = myMATH::linspace(_param.limit.x_min, _param.limit.x_max, _param.nx);
    _grid_y = myMATH::linspace(_param.limit.x_min, _param.limit.x_max, _param.nx);
    _grid_z = myMATH::linspace(_param.limit.x_min, _param.limit.x_max, _param.nx);

    // set parameters for displacement statistics
    _nBin_tot = _param.nBin_x * _param.nBin_y * _param.nBin_z;

    _m_xyz.resize(3);
    _c_xyz.resize(3);

    _m_xyz[0] = (_param.nBin_x-1)/(4*_param.r);
    _c_xyz[0] = (_param.nBin_x-1)/2;

    _m_xyz[1] = (_param.nBin_y-1)/(4*_param.r);
    _c_xyz[1] = (_param.nBin_y-1)/2;

    _m_xyz[2] = (_param.nBin_z-1)/(4*_param.r);
    _c_xyz[2] = (_param.nBin_z-1)/2;
}


template<class T3D>
void PredField::calDispField(std::vector<T3D> const& obj3d_list_prev, std::vector<T3D> const& obj3d_list_curr)
{
    _disp_field = Matrix<double>(_nGrid_tot,3,0);

    // calculate displacement field
    #pragma omp parallel for collapse(3)
    for (int x_id = 0; x_id < _param.nx; x_id ++)
    {
        for (int y_id = 0; y_id < _param.ny; y_id ++)
        {
            for (int z_id = 0; z_id < _param.nz; z_id ++)
            {
                std::vector<int> xyz_id = {x_id, y_id, z_id};
                calDispGrid(xyz_id, obj3d_list_prev, obj3d_list_curr);
            }
        }
    }
}


// For Pt3D list
void PredField::calDispField(std::vector<Pt3D> const& pt3d_list_prev, std::vector<Pt3D> const& pt3d_list_curr)
{
    _disp_field = Matrix<double>(_nGrid_tot,3,0);

    // calculate displacement field
    #pragma omp parallel for collapse(3)
    for (int x_id = 0; x_id < _param.nx; x_id ++)
    {
        for (int y_id = 0; y_id < _param.ny; y_id ++)
        {
            for (int z_id = 0; z_id < _param.nz; z_id ++)
            {
                std::vector<int> xyz_id = {x_id, y_id, z_id};
                calDispGrid(xyz_id, pt3d_list_prev, pt3d_list_curr);
            }
        }
    }

}


template<class T3D>
void PredField::calDispGrid(std::vector<int> const& xyz_id, std::vector<T3D> const& obj3d_list_prev, std::vector<T3D> const& obj3d_list_curr)
{
    int grid_id = mapGridID(xyz_id[0],xyz_id[1],xyz_id[2]);

    // find obj that is close to the grid
    std::vector<int> objID_list_prev, objID_list_curr;
    findNeighbor(objID_list_prev, xyz_id, obj3d_list_prev);
    findNeighbor(objID_list_curr, xyz_id, obj3d_list_curr);

    // if not find obj in one frame, then stop
    if (objID_list_curr.size()==0 || objID_list_curr.size()==0)
    {
        _disp_field(grid_id, 0) = 0;
        _disp_field(grid_id, 1) = 0;
        _disp_field(grid_id, 2) = 0;
        return;
    }

    // getting the displacement vectors (if there are particles in the search volume for both the frames)
    std::vector<std::vector<double>> disp_list; // all dx,dy,dz
    std::vector<double> disp(3,0); // dx,dy,dz
    int id_curr, id_prev;
    
    for (int j = 0; j < objID_list_curr.size(); j ++) 
    {
        id_curr = objID_list_curr[j];
        for (int k = 0; k < objID_list_prev.size(); k ++) 
        {  
            id_prev = objID_list_prev[k];

            disp[0] = obj3d_list_curr[id_curr]._pt_center[0] 
                    - obj3d_list_prev[id_prev]._pt_center[0];

            disp[1] = obj3d_list_curr[id_curr]._pt_center[1] 
                    - obj3d_list_prev[id_prev]._pt_center[1];

            disp[2] = obj3d_list_curr[id_curr]._pt_center[2] 
                    - obj3d_list_prev[id_prev]._pt_center[2];

            disp_list.push_back(disp);
        }
    }
    
    // update displacement pdf
    std::vector<double> disp_pdf(_nBin_tot);
    updateDispPDF(disp_pdf, disp_list);

    // find the peak location of the displacement pdf
    std::vector<double> disp_opt;
    findPDFPeakLoc(disp_opt, disp_pdf);

    _disp_field(grid_id, 0) = disp_opt[0];
    _disp_field(grid_id, 1) = disp_opt[1];
    _disp_field(grid_id, 2) = disp_opt[2];
}


// For Pt3D list
void PredField::calDispGrid(std::vector<int> const& xyz_id, std::vector<Pt3D> const& obj3d_list_prev, std::vector<Pt3D> const& obj3d_list_curr)
{
    int grid_id = mapGridID(xyz_id[0],xyz_id[1],xyz_id[2]);

    // find obj that is close to the grid
    std::vector<int> objID_list_prev, objID_list_curr;
    findNeighbor(objID_list_prev, xyz_id, obj3d_list_prev);
    findNeighbor(objID_list_curr, xyz_id, obj3d_list_curr);

    // if not find obj in one frame, then stop
    if (objID_list_curr.size()==0 || objID_list_curr.size()==0)
    {
        _disp_field(grid_id, 0) = 0;
        _disp_field(grid_id, 1) = 0;
        _disp_field(grid_id, 2) = 0;
        return;
    }

    // getting the displacement vectors (if there are particles in the search volume for both the frames)
    std::vector<std::vector<double>> disp_list; // all dx,dy,dz
    std::vector<double> disp(3,0); // dx,dy,dz
    int id_curr, id_prev;
    
    for (int j = 0; j < objID_list_curr.size(); j ++) 
    {
        id_curr = objID_list_curr[j];
        for (int k = 0; k < objID_list_prev.size(); k ++) 
        {  
            id_prev = objID_list_prev[k];

            disp[0] = obj3d_list_curr[id_curr][0] 
                    - obj3d_list_prev[id_prev][0];

            disp[1] = obj3d_list_curr[id_curr][1] 
                    - obj3d_list_prev[id_prev][1];

            disp[2] = obj3d_list_curr[id_curr][2] 
                    - obj3d_list_prev[id_prev][2];

            disp_list.push_back(disp);
        }
    }
    
    // update displacement pdf
    std::vector<double> disp_pdf(_nBin_tot);
    updateDispPDF(disp_pdf, disp_list);

    // find the peak location of the displacement pdf
    std::vector<double> disp_opt;
    findPDFPeakLoc(disp_opt, disp_pdf);

    _disp_field(grid_id, 0) = disp_opt[0];
    _disp_field(grid_id, 1) = disp_opt[1];
    _disp_field(grid_id, 2) = disp_opt[2];
}


template<class T3D>
void PredField::findNeighbor (std::vector<int>& pt_list_id, std::vector<int> const& xyz_id, std::vector<T3D> const& obj3d_list)
{
    double x = _grid_x[xyz_id[0]];
    double y = _grid_y[xyz_id[1]];
    double z = _grid_z[xyz_id[2]];

    double dx, dy, dz, distsqr;
    double rsqr = _param.r*_param.r;

    for (int i = 0; i < obj3d_list.size(); i++)
    {
        dx = obj3d_list[i]._pt_center[0] - x;
        dy = obj3d_list[i]._pt_center[1] - y;
        dz = obj3d_list[i]._pt_center[2] - z;

        if (dx<_param.r && dx>-_param.r &&
            dy<_param.r && dy>-_param.r &&
            dz<_param.r && dz>-_param.r)
        {
            distsqr = dx*dx + dy*dy + dz*dz;
            if (distsqr < rsqr)
            {
                pt_list_id.push_back(i);
            }
        }
    }
}


void PredField::findNeighbor (std::vector<int>& pt_list_id, std::vector<int> const& xyz_id, std::vector<Pt3D> const& pt3d_list)
{
    double x = _grid_x[xyz_id[0]];
    double y = _grid_y[xyz_id[1]];
    double z = _grid_z[xyz_id[2]];

    double dx, dy, dz, distsqr;
    double rsqr = _param.r*_param.r;

    for (int i = 0; i < pt3d_list.size(); i++)
    {
        dx = pt3d_list[i][0] - x;
        dy = pt3d_list[i][1] - y;
        dz = pt3d_list[i][2] - z;

        if (dx<_param.r && dx>-_param.r &&
            dy<_param.r && dy>-_param.r &&
            dz<_param.r && dz>-_param.r)
        {
            distsqr = dx*dx + dy*dy + dz*dz;
            if (distsqr < rsqr)
            {
                pt_list_id.push_back(i);
            }
        }
    }
}


int PredField::mapGridID (int x_id, int y_id, int z_id)
{
    return x_id * _param.ny*_param.nz + y_id * _param.nz + z_id;
}


int PredField::mapBinID (int xBin_id, int yBin_id, int zBin_id)
{
    return xBin_id*_param.nBin_y*_param.nBin_z + yBin_id*_param.nBin_z + zBin_id;
}


void PredField::updateDispPDF(std::vector<double>& disp_pdf, std::vector<std::vector<double>> const& disp_list)
{
    if (disp_pdf.size() != _nBin_tot)
    {
        disp_pdf.resize(_nBin_tot);
    }

    // initialize disp_pdf
    std::fill(disp_pdf.begin(), disp_pdf.end(), 0);

    for (int id = 0; id < disp_list.size(); id ++)
    {
        double xBin = _m_xyz[0]*disp_list[id][0] + _c_xyz[0];
        double yBin = _m_xyz[1]*disp_list[id][1] + _c_xyz[1];
        double zBin = _m_xyz[2]*disp_list[id][2] + _c_xyz[2];

        int xBin_id = std::round(xBin);
        int yBin_id = std::round(yBin);
        int zBin_id = std::round(zBin);

        for (int i = std::max(0, xBin_id-1); 
             i <= std::min(_param.nBin_x-1, xBin_id+1); 
             i ++) 
        {
            double xx = i - xBin;
            for (int j = std::max(0, yBin_id-1); 
                 j <= std::min(_param.nBin_y-1, yBin_id+1); 
                 j ++) 
            {
                double yy = j - yBin;
                for (int k = std::max(0, zBin_id-1); 
                     k <= std::min(_param.nBin_z-1, zBin_id+1); 
                     k ++) 
                {
                    double zz = k - zBin;
                    int bin_id = mapBinID(i,j,k);
                    disp_pdf[bin_id] += std::exp(-(xx*xx + yy*yy + zz*zz));
                }
            }
        }
    }
}


void PredField::findPDFPeakLoc(std::vector<double>& dist_opt, std::vector<double> const& disp_pdf)
{
    if (dist_opt.size() != 3)
    {
        dist_opt.resize(3);
    }

    std::vector<int> peak_id;
    findPDFPeakID(peak_id, disp_pdf);

    // finding the peak using 1D Gaussian in x, y & z direction
    int x1, y1, z1, x2, y2, z2, x3, y3, z3;
    for (int i = 0; i < 3; i ++)
    {
        if (peak_id[i] == 0)
        {
            dist_opt[i] = 0.0; // no need to calculate, it is possible that there is only one bin 
            // dist_opt[i] = -2*_param.r;
        }
        else if (peak_id[i] == _param.nBin_x-1)
        {
            dist_opt[i] = 2*_param.r;
        }
        else
        {
            x1 = peak_id[i] - 1 * (i==0);
            y1 = peak_id[i] - 1 * (i==1);
            z1 = peak_id[i] - 1 * (i==2);

            x2 = peak_id[i];
            y2 = peak_id[i];
            z2 = peak_id[i];

            x3 = peak_id[i] + 1 * (i==0);
            y3 = peak_id[i] + 1 * (i==1);
            z3 = peak_id[i] + 1 * (i==2);

            double loc = fitPeakBinLocGauss(
                peak_id[i]-1, disp_pdf[mapBinID(x1,y1,z1)],
                peak_id[i],   disp_pdf[mapBinID(x2,y2,z2)],
                peak_id[i]+1, disp_pdf[mapBinID(x3,y3,z3)]
            );

            dist_opt[i] = (loc - _c_xyz[i]) / _m_xyz[i];
        }
    }
}


void PredField::findPDFPeakID(std::vector<int>& peak_id, std::vector<double> const& disp_pdf)
{
    if (peak_id.size() != 3)
    {
        peak_id.resize(3);
    }

    std::fill(peak_id.begin(), peak_id.end(), 0);

    double val = disp_pdf[0];
    int bin_id = 0;
    for (int i = 0; i < _param.nBin_x; i ++)
    {
        for (int j = 0; j < _param.nBin_y; j ++)
        {
            for (int k = 0; k < _param.nBin_z; k ++)
            {
                bin_id = mapBinID(i,j,k);
                if (disp_pdf[bin_id] > val)
                {
                    val = disp_pdf[bin_id];

                    peak_id[0] = i;
                    peak_id[1] = j;
                    peak_id[2] = k;
                }
            }
        }
    }
}


double PredField::fitPeakBinLocGauss (double y1, double v1, double y2, double v2, double y3, double v3)
{
    double ln_z1, ln_z2, ln_z3;

    ln_z1 = v1 < LOGSMALLNUMBER ? std::log(LOGSMALLNUMBER) : std::log(v1);
    ln_z2 = v2 < LOGSMALLNUMBER ? std::log(LOGSMALLNUMBER) : std::log(v2);
    ln_z3 = v3 < LOGSMALLNUMBER ? std::log(LOGSMALLNUMBER) : std::log(v3);

    double yc = -0.5 * (  (ln_z1 * double((y2 * y2) - (y3 * y3))) 
                        - (ln_z2 * double((y1 * y1) - (y3 * y3))) 
                        + (ln_z3 * double((y1 * y1) - (y2 * y2))) ) 
                     / (  (ln_z1 * double(y3 - y2)) 
                        - (ln_z3 * double(y1 - y2)) 
                        + (ln_z2 * double(y1 - y3)) );

    if (!std::isfinite(yc))
    {
        yc = y2;
    }

    return yc;
}


#endif