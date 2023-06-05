#include "PredField.h"

void PredField::SetGrid ()
{
    int nx = _n_xyz[0];
    int ny = _n_xyz[1];
    int nz = _n_xyz[2];

    std::vector<double> grid_x = myMATH::Linspace(_limit._x_min, _limit._x_max, nx);
    std::vector<double> grid_y = myMATH::Linspace(_limit._y_min, _limit._y_max, ny);
    std::vector<double> grid_z = myMATH::Linspace(_limit._z_min, _limit._z_max, nz);

    int id = 0;
    for (int i = 0; i < nx; i ++)
    {
        for (int j = 0; j < ny; j ++)
        {
            for (int k = 0; k < nz; k ++)
            {
                _grid(0, id) = grid_x[i];
                _grid(1, id) = grid_y[j];
                _grid(2, id) = grid_z[k];
                id ++;
            }
        }
    }
}


void PredField::Field()
{
    double rsqr = pow(_r, 2);


    #pragma omp parallel //num_threads(8)
    {
        #pragma omp for
        for (int i = 0; i < _n_tot; i++) 
        {
            std::vector<int> prev_id;
            std::vector<int> curr_id;
            std::vector<std::vector<double>> displacements;
            
            // getting the points within the interrogation sphere around the grid point
            // previous frame
            FindVolPt(prev_id, rsqr, i, PREV_FRAME);
            // current frame
            FindVolPt(curr_id, rsqr, i, CURR_FRAME);

            // getting the displacement vectors (if there are particles in the search volume for both the frames)
            std::vector<double> disp(4, 0); // dx,dy,dz,I_curr*I_prev
            int curr, prev;
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
                std::vector<std::vector<std::vector<int>>> disp_map(3, std::vector<std::vector<int>> (3, std::vector<int> (3,0)));

                // getting the 3D displacement correlation map
                DispMap(disp_map, displacements);

                // finding the peak of displacement map
                std::vector<double> peak(3);
                peak = DispMapPeak(disp_map);
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

}
