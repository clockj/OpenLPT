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
void PredField<T>::SetPtList ()
{
    for (int i = 0; i < _obj_list_prev.size(); i ++)
    {
        _pt_list_prev[i] = _obj_list_prev[i].GetCenterPos();
    }

    for (int i = 0; i < _obj_list_curr.size(); i ++)
    {
        _pt_list_curr[i] = _obj_list_curr[i].GetCenterPos();
    }
}



template<class T>
void PredField<T>::Field()
{
    double rsqr = _r*_r;

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

            // // check id=[3500,5157,6291,6809,8446]-1
            // if (i==3499 || i==5156 || i==6290 || i==6808 || i==8445)
            // {
            //     std::cout << i << ";" << _grid(0,i) << "," << _grid(1,i) << "," << _grid(2,i) << std::endl;
            //     std::ofstream file_prev_id("prev_pt3d_"+std::to_string(i)+".csv", std::ios::out);
            //     if (prev_id.size()!=0)
            //     {
            //         file_prev_id.precision(8);
            //         for (int j = 0; j < prev_id.size(); j ++)
            //         {
            //             file_prev_id << prev_id[j] << ",";
            //             file_prev_id << _pt_list_prev[prev_id[j]](0,0) << "," 
            //                          << _pt_list_prev[prev_id[j]](1,0) << ","
            //                          << _pt_list_prev[prev_id[j]](2,0) << "\n";
            //         }
            //     }
            //     else 
            //     {
            //         file_prev_id << "\n";
            //     }
            //     file_prev_id.close();

            //     std::ofstream file_curr_id("curr_pt3d_"+std::to_string(i)+".csv", std::ios::out);
            //     file_curr_id.precision(8);
            //     if (curr_id.size()!=0)
            //     {
            //         for (int j = 0; j < curr_id.size(); j ++)
            //         {
            //             file_curr_id << curr_id[j] << ",";
            //             file_curr_id << _pt_list_curr[curr_id[j]](0,0) << "," 
            //                          << _pt_list_curr[curr_id[j]](1,0) << ","
            //                          << _pt_list_curr[curr_id[j]](2,0) << "\n";
            //         }
            //     }
            //     else
            //     {
            //         file_curr_id << "\n";
            //     }
            //     file_curr_id.close();
            // }

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

                        disp[0] = _pt_list_curr[curr](0,0) - _pt_list_prev[prev](0,0);
                        disp[1] = _pt_list_curr[curr](1,0) - _pt_list_prev[prev](1,0);
                        disp[2] = _pt_list_curr[curr](2,0) - _pt_list_prev[prev](2,0);
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

    // // rule out outlier
    // #pragma omp parallel
    // {
    //     #pragma omp for
    //     for (int i = 0; i < 3; i ++)
    //     {
    //         std::vector<double> row;
    //         std::vector<int> grid_id;
    //         for (int j = 0; j < _n_tot; j ++)
    //         {
    //             if (_disp_field(i,j)!=0)
    //             {
    //                 row.push_back(_disp_field(i,j));
    //                 grid_id.push_back(j);
    //             }
    //         }

    //         if  (row.size() == 0)
    //         {
    //             continue;
    //         }
    //         std::vector<int> is_outlier(row.size(),0);
    //         myMATH::IsOutlier(is_outlier, row);

    //         for (int j = 0; j < is_outlier.size(); j ++)
    //         {
    //             if (is_outlier[j]) // if it is an outlier
    //             {
    //                 _disp_field(0, grid_id[j]) = 0;
    //                 _disp_field(1, grid_id[j]) = 0;
    //                 _disp_field(2, grid_id[j]) = 0;
    //             }
    //         }
    //     }
    // }

    t_end = clock();
    std::cout << "Displacement field setup time: " 
              << (double) (t_end - t_start)/CLOCKS_PER_SEC;

}


template<class T>
void PredField<T>::FindVolPt(std::deque<int>& pt_list_id, double rsqr, int grid_id, FrameTypeID frame)
{
    double x = _grid(0, grid_id);
    double y = _grid(1, grid_id);
    double z = _grid(2, grid_id);
    double dx, dy, dz, distsqr;
    
    // bool is_save = false;
    // std::ofstream file;
    if (frame == PREV_FRAME)
    {
        // if (grid_id==3499 || grid_id==5156 || grid_id==6290 || grid_id==6808 || grid_id==8445)
        // {
        //     std::cout << "FindVolPt,prev: " << grid_id << ";" << _grid(0,grid_id) << "," << _grid(1,grid_id) << "," << _grid(2,grid_id) << std::endl;
        //     is_save = true;
        //     file.open("prev_search_"+std::to_string(grid_id)+".csv", std::ios::out);
        //     file.precision(8);
        // }

        for (int i = 0; i < _pt_list_prev.size(); i ++)
        {   
            dx = _pt_list_prev[i](0,0) - x;
            dy = _pt_list_prev[i](1,0) - y;
            dz = _pt_list_prev[i](2,0) - z;
            
            if (dx<_r && dx>-_r &&
                dy<_r && dy>-_r &&
                dz<_r && dz>-_r)
            {
                distsqr = dx*dx + dy*dy + dz*dz;
                if (distsqr < rsqr)
                {
                    pt_list_id.push_back(i);
                }
            }

            // if (is_save)
            // {
            //     file << i << "," << _pt_list_prev[i](0,0) << "," 
            //                      << _pt_list_prev[i](1,0) << "," 
            //                      << _pt_list_prev[i](2,0) << ","
            //                      << dx << "," << dy << "," << dz << "," 
            //                      << distsqr << "," << rsqr << "\n";
            // }
        }
        // if (is_save)
        // {
        //     file.close();
        // }
    }
    else if (frame == CURR_FRAME)
    {
        // if (grid_id==3499 || grid_id==5156 || grid_id==6290 || grid_id==6808 || grid_id==8445)
        // {
        //     std::cout << "FindVolPt,curr: " << grid_id << ";" << _grid(0,grid_id) << "," << _grid(1,grid_id) << "," << _grid(2,grid_id) << std::endl;
        //     is_save = true;
        //     file.open("curr_search_"+std::to_string(grid_id)+".csv", std::ios::out);
        //     file.precision(8);
        // }

        for (int i = 0; i < _pt_list_curr.size(); i ++)
        {   
            dx = _pt_list_curr[i](0,0) - x;
            dy = _pt_list_curr[i](1,0) - y;
            dz = _pt_list_curr[i](2,0) - z;
            if (dx<_r && dx>-_r &&
                dy<_r && dy>-_r &&
                dz<_r && dz>-_r)
            {
                distsqr = dx*dx + dy*dy + dz*dz;
                if (distsqr < rsqr)
                {
                    pt_list_id.push_back(i);
                }
            }

            // if (is_save)
            // {
            //     file << i << "," << _pt_list_curr[i](0,0) << "," 
            //                      << _pt_list_curr[i](1,0) << "," 
            //                      << _pt_list_curr[i](2,0) << ","
            //                      << dx << "," << dy << "," << dz << "," 
            //                      << distsqr << "," << rsqr << "\n";
            // }
        }
        // if (is_save)
        // {
        //     file.close();
        // }
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
        peak[0] = _c;
        // peak[0] = 0;
    }
    else if (x == _size - 1)
    {
        peak[0] = 2*_r;
        // peak[0] = _size-1;
    }   
    else
    {
        peak[0] = Gauss1DPeak(x - 1, disp_map[x - 1][y][z], x, disp_map[x][y][z], x + 1, disp_map[x + 1][y][z]);
    }

    if (y == 0)
    {
        peak[1] = _c;
        // peak[1] = 0;
    }
    else if (y == _size - 1)
    {
        peak[1] = 2*_r;
        // peak[1] = _size - 1;
    }
    else
    {
        peak[1] = Gauss1DPeak(y - 1, disp_map[x][y - 1][z], y, disp_map[x][y][z], y + 1, disp_map[x][y + 1][z]);
    }

    if (z == 0)
    {
        peak[2] = _c;
        // peak[2] = 0;
    }
    else if (z == _size - 1)
    {
        peak[2] = 2*_r;
        // peak[2] = _size - 1;
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
        lnz1 = std::log(0.0001);
    }
    else 
    {
        lnz1 = std::log(v1);
    }

    if (v2 == 0) 
    {
        lnz2 = std::log(0.0001);
    }
    else 
    {
        lnz2 = std::log(v2);
    }

    if (v3 == 0) 
    {
        lnz3 = std::log(0.0001);
    }
    else 
    {
        lnz3 = std::log(v3);
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

    index_x = std::min(index_x, _n_xyz[0]-2);
    index_y = std::min(index_y, _n_xyz[1]-2);
    index_z = std::min(index_z, _n_xyz[2]-2);

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