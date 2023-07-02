#ifndef PREDFIELD_HPP
#define PREDFIELD_HPP

#include "PredField.h"

template<class T>
void PredField<T>::SetGrid ()
{
    int nx = _n_xyz[0];
    int ny = _n_xyz[1];
    int nz = _n_xyz[2];

    _dx = (_limit._x_max-_limit._x_min) / (nx-1);
    _dy = (_limit._y_max-_limit._y_min) / (ny-1);
    _dz = (_limit._z_max-_limit._z_min) / (nz-1);

    _grid_x = myMATH::Linspace(_limit._x_min, _limit._x_max, nx);
    _grid_y = myMATH::Linspace(_limit._y_min, _limit._y_max, ny);
    _grid_z = myMATH::Linspace(_limit._z_min, _limit._z_max, nz);

    int id = 0;
    for (int i = 0; i < nx; i ++)
    {
        for (int j = 0; j < ny; j ++)
        {
            for (int k = 0; k < nz; k ++)
            {
                _grid(0, id) = _grid_x[i];
                _grid(1, id) = _grid_y[j];
                _grid(2, id) = _grid_z[k];
                id ++;
            }
        }
    }
}


template<class T>
void PredField<T>::Field()
{
    double rsqr = pow(_r, 2);

    clock_t t_start, t_end;
    t_start = clock();

    // omp_set_num_threads(1);
    #pragma omp parallel //num_threads(8)
    {
        #pragma omp for
        for (int i = 0; i < _n_tot; i ++) 
        {
            std::deque<int> prev_id;
            std::deque<int> curr_id;
            std::vector<std::vector<double>> displacements;
            
            // getting the points within the interrogation sphere around the grid point
            // previous frame
            FindVolPt(prev_id, rsqr, i, PREV_FRAME);
            // current frame
            FindVolPt(curr_id, rsqr, i, CURR_FRAME);

            // getting the displacement vectors (if there are particles in the search volume for both the frames)
            std::vector<double> disp(4, 0); // dx,dy,dz,I_curr*I_prev
            int curr, prev;
            Matrix<double> pt_pre(3,1,0);
            Matrix<double> pt_cur(3,1,0);
            if (curr_id.size() != 0 && prev_id.size() != 0) 
            {
                for (int j = 0; j < curr_id.size(); j ++) 
                {
                    curr = curr_id[j];
                    for (int k = 0; k < prev_id.size(); k ++) 
                    {  
                        prev = prev_id[k];
                        pt_pre = _pt_list_prev[prev].GetCenterPos();
                        pt_cur = _pt_list_curr[curr].GetCenterPos();

                        disp[0] = pt_cur(0,0) - pt_pre(0,0);
                        disp[1] = pt_cur(1,0) - pt_pre(1,0);
                        disp[2] = pt_cur(2,0) - pt_pre(2,0);
                        disp[3] = 1; // currFrame[curr]->Info() * prevFrame[prev]->Info();
                        displacements.push_back(disp);
                    }
                }

                // intializing the displacement map to 0
                std::vector<std::vector<std::vector<double>>> disp_map(_size, std::vector<std::vector<double>> (_size, std::vector<double> (_size, 0)));

                // getting the 3D displacement correlation map
                DispMap(disp_map, displacements);

                // finding the peak of displacement map
                std::vector<double> peak = DispMapPeak(disp_map);

                // saving the peak dx,dy,dz to field
                _disp_field(0, i) = peak[0];
                _disp_field(1, i) = peak[1];
                _disp_field(2, i) = peak[2];
            }
            else 
            {
                _disp_field(0, i) = 0;
                _disp_field(1, i) = 0;
                _disp_field(2, i) = 0;
            }
        }
    }

    t_end = clock();
    std::cout << "Displacement field setup time: " 
              << (double) (t_end - t_start)/CLOCKS_PER_SEC
              << std::endl;

}


template<class T>
void PredField<T>::FindVolPt(std::deque<int>& pt_list_id, double rsqr, int grid_id, FrameTypeID frame)
{
    double x = _grid(0, grid_id);
    double y = _grid(1, grid_id);
    double z = _grid(2, grid_id);
    double x_temp, y_temp, z_temp, distsqr;
    
    if (frame == PREV_FRAME)
    {
        for (int i = 0; i < _pt_list_prev.size(); i ++)
        {   
            x_temp = _pt_list_prev[i].GetCenterPos()(0,0);
            y_temp = _pt_list_prev[i].GetCenterPos()(1,0);
            z_temp = _pt_list_prev[i].GetCenterPos()(2,0);
            if (x_temp<x+_r && x_temp>x-_r &&
                y_temp<y+_r && y_temp>y-_r &&
                z_temp<z+_r && z_temp>z-_r)
            {
                distsqr = std::pow(x_temp-x,2) + std::pow(y_temp-y,2) + std::pow(z_temp-z,2);
                if (distsqr < rsqr)
                {
                    pt_list_id.push_back(i);
                }
            }
        }
    }
    else if (frame == CURR_FRAME)
    {
        for (int i = 0; i < _pt_list_curr.size(); i ++)
        {   
            x_temp = _pt_list_curr[i].GetCenterPos()(0,0);
            y_temp = _pt_list_curr[i].GetCenterPos()(1,0);
            z_temp = _pt_list_curr[i].GetCenterPos()(2,0);
            if (x_temp<x+_r && x_temp>x-_r &&
                y_temp<y+_r && y_temp>y-_r &&
                z_temp<z+_r && z_temp>z-_r)
            {
                distsqr = std::pow(x_temp-x,2) + std::pow(y_temp-y,2) + std::pow(z_temp-z,2);
                if (distsqr < rsqr)
                {
                    pt_list_id.push_back(i);
                }
            }
        }
    }
    else 
    {
        std::cerr << "PredField::FindVolPt: frame type wrong!" << std::endl;
        throw error_type;
    }
}


template <class T>
void PredField<T>::DispMap(std::vector<std::vector<std::vector<double>>>& disp_map, std::vector<std::vector<double>> const& displacements)
{
    for (int d = 0; d < displacements.size(); d ++) 
    {
        double dx = _m * displacements[d][0] + _c;
        double dy = _m * displacements[d][1] + _c; 
        double dz = _m * displacements[d][2] + _c;
        double I = displacements[d][3];
        int x = round(dx), y = round(dy), z = round(dz);
        for (int i = std::max(0, x-1); i <= std::min(_size-1, x+1); i ++) 
        {
            double xx = i - dx;
            for (int j = std::max(0, y-1); j <= std::min(_size-1, y+1); j ++) 
            {
                double yy = j - dy;
                for (int k = std::max(0, z-1); k <= std::min(_size-1, z+1); k ++) {
                    double zz = k - dz;
                    disp_map[i][j][k] += I*std::exp(-(xx*xx + yy*yy + zz*zz));
                }
            }
        }	
    }
}


template <class T>
std::vector<double> PredField<T>::DispMapPeak(std::vector<std::vector<std::vector<double>>>& disp_map)
{
    // finding the index of largest element
    std::vector<int> index(3,0);
    index = MaxElemID(disp_map);
    int x = index[0];
    int y = index[1];
    int z = index[2];

    // finding the peak using 1D Gaussian in x, y & z direction
    std::vector<double> peak(3);
    if (x == 0)
    {
        peak[0] = 0;
    }
    else if (x == _size - 1)
    {
        peak[0] = _size - 1;
    }   
    else
    {
        peak[0] = Gauss1DPeak(x - 1, disp_map[x - 1][y][z], x, disp_map[x][y][z], x + 1, disp_map[x + 1][y][z]);
    }

    if (y == 0)
    {
        peak[1] = 0;
    }
    else if (y == _size - 1)
    {
        peak[1] = _size - 1;
    }
    else
    {
        peak[1] = Gauss1DPeak(y - 1, disp_map[x][y - 1][z], y, disp_map[x][y][z], y + 1, disp_map[x][y + 1][z]);
    }

    if (z == 0)
    {
        peak[2] = 0;
    }
    else if (z == _size - 1)
    {
        peak[2] = _size - 1;
    }
    else
    {
        peak[2] = Gauss1DPeak(z - 1, disp_map[x][y][z - 1], z, disp_map[x][y][z], z + 1, disp_map[x][y][z + 1]);
    }

    return peak;
}


template <class T>
std::vector<int> PredField<T>::MaxElemID(std::vector<std::vector<std::vector<double>>>& disp_map)
{
    std::vector<int> max_id (3,0);

    double val = disp_map[0][0][0];
    for (int i = 0; i < _size; i ++)
    {
        for (int j = 0; j < _size; j ++)
        {
            for (int k = 0; k < _size; k ++)
            {
                if (disp_map[i][j][k] > val)
                {
                    val = disp_map[i][j][k];
                    max_id[0] = i; 
                    max_id[1] = j;
                    max_id[2] = k;
                }
            }
        }
    }

    return max_id;
}
    

template <class T>
double PredField<T>::Gauss1DPeak(double y1, double v1, double y2, double v2, double y3, double v3)
{
    double lnz1, lnz2, lnz3;

    if (v1 == 0) 
    {
        lnz1 = log(0.0001);
    }
    else 
    {
        lnz1 = log(v1);
    }

    if (v2 == 0) 
    {
        lnz2 = log(0.0001);
    }
    else 
    {
        lnz2 = log(v2);
    }

    if (v3 == 0) 
    {
        lnz3 = log(0.0001);
    }
    else 
    {
        lnz3 = log(v3);
    }
    
    double yc = - 0.5 * ( lnz1*(y2*y2-y3*y3) - lnz2*(y1*y1-y3*y3) + lnz3*(y1*y1-y2*y2) )
                      / ( lnz1*(y3-y2) - lnz3*(y1-y2) + lnz2*(y1-y3) );

    //cout << "xc: " << xc << " yc: " << yc << endl;
    return (yc - _c)/_m;
}


template <class T>
Matrix<double> PredField<T>::PtInterp(Matrix<double> const& pt)
{
    // find out the limits of interpolation cube 
    double pt_x = std::min(
        std::max(pt(0,0), _limit._x_min),
        _limit._x_max
    );
    double pt_y = std::min(
        std::max(pt(1,0), _limit._y_min),
        _limit._y_max
    );
    double pt_z = std::min(
        std::max(pt(2,0), _limit._z_min),
        _limit._z_max
    );

    int index_x = std::max(
        0, 
        (int) std::floor((pt_x - _limit._x_min) / _dx)
    );
    int index_y = std::max(
        0, 
        (int) std::floor((pt_y - _limit._y_min) / _dy)
    );
    int index_z = std::max(
        0, 
        (int) std::floor((pt_z - _limit._z_min) / _dz)
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

    Matrix<double> disp(3,1,0);
    std::vector<double> pt_vec = {pt_x, pt_y, pt_z};

    for (int j = 0; j < 3; j ++)
    {
        std::vector<double> field = {
            _disp_field(j, i_000),
            _disp_field(j, i_100),
            _disp_field(j, i_101),
            _disp_field(j, i_001),
            _disp_field(j, i_010),
            _disp_field(j, i_110),
            _disp_field(j, i_111),
            _disp_field(j, i_011)
        };
        disp(j,0) = myMATH::TriLinearInterp(grid_limit, field, pt_vec);
    }

    return disp;
}



#endif