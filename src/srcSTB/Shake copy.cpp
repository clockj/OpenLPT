#ifndef SHAKE_HPP
#define SHAKE_HPP

#include "Shake.h"

PixelRange Shake::FindSearchRegion (int cam_use_id, int center_row, int center_col, int half_width_pixel)
{
    // int center_row = tracer.tracer_center_pixel(1,0);
    // int center_col = tracer.tracer_center_pixel(0,0);
    int n_row = _cam_use[cam_use_id].GetNpixh();
    int n_col = _cam_use[cam_use_id].GetNpixw();

    PixelRange search_range;
    search_range._row_min = std::min(
        std::max(
            0,
            center_row - half_width_pixel
        ), 
    n_row-1);
    search_range._col_min = std::min(
        std::max(
            0,
            center_col - half_width_pixel
        ),
    n_col-1);
    search_range._row_max = std::max(
        std::min(
            n_row,
            center_row + half_width_pixel + 1
        ),
    1);
    search_range._col_max = std::max(
        std::min(
            n_col,
            center_col + half_width_pixel + 1
        ),
    1);

    return search_range;
}


// (x, y) = (col_id, row_id)
double Shake::GaussianProjection (Matrix<double>& tracer_center_pixel, std::vector<double>& otf_param, int x, int y)
{
    double xx =   (double(x) - tracer_center_pixel(0,0)) * std::cos(otf_param[3]) 
                + (double(y) - tracer_center_pixel(1,0)) * std::sin(otf_param[3]);
    double yy = - (double(x) - tracer_center_pixel(0,0)) * std::sin(otf_param[3]) 
                + (double(y) - tracer_center_pixel(1,0)) * std::cos(otf_param[3]);
    double value = otf_param[0] 
                   * std::exp( - (otf_param[1] * xx*xx + otf_param[2] * yy*yy) );
    return std::fabs(value);
}


double Shake::CalPointResidue (Matrix<double>& pos_new, std::vector<PixelRange>& search_range_list, std::vector<Matrix<double>>& aug_img_list)
{
    double residue = 0;

    Matrix<double> pos_2d_pixel(3,1,0);
    int id;
    for (int m = 0; m < _n_cam_use; m ++)
    {
        pos_2d_pixel = _cam_use[m].ImgMMToPixel(
            _cam_use[m].WorldToImgMM(pos_new)
        );

        id = _cam_id[m];
        std::vector<double> otf_param = _otf.GetOTFParam(id, pos_new);

        // Calculating the residual for updated 3D position
        // in the original code, 
        //  residue = sum (sum residue^2) for all cam 
        // should it be 
        //  residue = sum sqrt(sum residue^2) / _n_cam_use 
        PixelRange search_range = search_range_list[m];
        int i = 0; 
        for (int row_id = search_range._row_min; row_id < search_range._row_max; row_id ++)
        {
            int j = 0; 
            for (int col_id = search_range._col_min; col_id < search_range._col_max; col_id ++)
            {
                residue += std::pow(
                    aug_img_list[m](i,j)
                    - GaussianProjection(
                            pos_2d_pixel,
                            otf_param,
                            col_id, row_id
                        ),
                    2
                );
                j ++;
            }
            i ++;
        }
    }

    // residue /= _n_cam_use;
    return residue;
}


// only use the initial 3 elements
void Shake::FitQuadratic (std::vector<double>& coeff, double* array, std::vector<double>& residue)
{
    // std::vector<double> coeff(3,0);

    // residue = coeff[0]*x^2 + coeff[1]*x + coeff[2]
    coeff[0] = (
                (residue[1]-residue[0]) / (array[1]-array[0]) 
                - (residue[2]-residue[0]) / (array[2]-array[0])
               ) / (array[1]-array[2]);

    coeff[1] = (residue[1]-residue[0]) / (array[1]-array[0]) 
                - coeff[0] * (array[1]+array[0]);

    coeff[2] = residue[2] - coeff[0]*array[2]*array[2] - coeff[1]*array[2];

    // return coeff;
}


Matrix<double> Shake::UpdatePos3D (Matrix<double> const& pos_old, std::vector<PixelRange>& search_range_list, std::vector<Matrix<double>>& aug_img_list, double delta, int half_width)
{
    Matrix<double> pos_new(pos_old);

    double delta_list[3];
    delta_list[0] = - delta;
    delta_list[1] = 0;
    delta_list[2] = delta;

    double array_list[4] = {0};
    std::vector<double> coeff(3, 0);
    std::vector<double> residue_list(4, 0);

    // shaking on x,y,z direction
    for (int i = 0; i < 3; i ++)
    {
        for (int j = 0; j < 3; j ++)
        {
            array_list[j] = pos_old(i,0) + delta_list[j];
            pos_new(i,0) = array_list[j];
            residue_list[j] = CalPointResidue(
                pos_new,
                search_range_list,
                aug_img_list
            );
        }
        FitQuadratic(coeff, array_list, residue_list);

        if (coeff[0] == 0)
        {
            if (coeff[1] == 0)
            {
                pos_new(i,0) = array_list[1];
            }
            else if (coeff[1] > 0)
            {
                pos_new(i,0) = array_list[0];
            }
            else
            {
                pos_new(i,0) = array_list[2];
            }
        }
        else 
        {
            array_list[3] = - coeff[1] / (2 * coeff[0]);
            if (array_list[3]>array_list[0] && array_list[3]<array_list[2])
            {
                pos_new(i,0) = array_list[3];
                residue_list[3] = CalPointResidue(
                    pos_new,
                    search_range_list,
                    aug_img_list
                );
            }
            else
            {
                residue_list[3] = residue_list[0] + residue_list[1] + residue_list[2] + 1; // set a maximum residue value
            }

            pos_new(i,0) = array_list[myMATH::MinID<double>(residue_list)];
        }
        
        // TODO: update Augemented image according to old position but new range?
        for (int m = 0; m < _n_cam_use; m ++)
        {
            Matrix<double> pos_2d_pixel(3,1,0);
            int id;
            for (int m = 0; m < _n_cam_use; m ++)
            {
                pos_2d_pixel = _cam_use[m].ImgMMToPixel(
                    _cam_use[m].WorldToImgMM(pos_new)
                );
                PixelRange search_range = FindSearchRegion(m, pos_2d_pixel(1,0), pos_2d_pixel(0,0), half_width);

                Matrix<double> aug_img(aug_img_list[m]);
                int n_col_new = search_range.GetNumOfCol();
                int n_row_new = search_range.GetNumOfRow();
                int n_col_old = search_range_list[m].GetNumOfCol();
                int n_row_old = search_range_list[m].GetNumOfRow();
                if (n_col_new < n_col_old || 
                    n_row_new < n_row_old)
                {
                    continue;
                }
                else if (n_col_new != n_col_old || 
                         n_row_new != n_row_old)
                {
                    aug_img_list[m] = Matrix<double> (n_row_new, n_col_new, 0);
                }
                
                // Get otf param
                id = _cam_id[m];
                std::vector<double> otf_param = _otf.GetOTFParam(id, pos_new);

                // update aug_img based on new range
                int row_min_old = search_range_list[m]._row_min;
                int row_max_old = search_range_list[m]._row_max;
                int col_min_old = search_range_list[m]._col_min;
                int col_max_old = search_range_list[m]._col_max;

                int i = 0; 
                for (int row_id = search_range._row_min; row_id < search_range._row_max; row_id ++)
                {
                    int j = 0; 
                    for (int col_id = search_range._col_min; col_id < search_range._col_max; col_id ++)
                    {
                        bool judge = row_id <  row_max_old && 
                                     col_id <  col_max_old &&
                                     row_id >= row_min_old &&
                                     col_id >= col_min_old;
                        if (judge)
                        {
                            aug_img_list[m](i,j) = aug_img(row_id-row_min_old, col_id-col_min_old);
                        }
                        else
                        {
                            aug_img_list[m](i,j) = _res_img_list[id](row_id, col_id);
                        }
                        j ++;
                    }
                    i ++;
                }

                // update search range
                search_range_list[m] = search_range;
            }
        }
    }

    return pos_new;
}


// remove search_range_small, directly use new range 
// if use new range, if it is out of the original range, use the original image intensity
double Shake<T>::CalObjectIntensity (TracerInfo& tracer, std::vector<PixelRange>& search_range_list, std::vector<Matrix<double>>& aug_img_list, double intensity, double min_ratio_of_median)
{
    Matrix<double> pos_new(tracer.GetCenterPos());

    // Pixel range on the cameras around the updated particle centers
    Matrix<double> tracer_center_pixel(3,1,0);
    std::vector<PixelRange> search_range_list_new(_n_cam_use);

    int n_outrange = 0;
    for (int m = 0; m < _n_cam_use; m ++)
    {
        tracer_center_pixel = _cam_use[m].ImgMMToPixel(
            _cam_use[m].WorldToImgMM(
                pos_new
            )
        );

        int row_id = tracer_center_pixel(1,0);
        int col_id = tracer_center_pixel(0,0);
        int n_row = _cam_use[m].GetHpix();
        int n_col = _cam_use[m].GetWpix();

        search_range_list_new[m] = FindSearchRegion(
            m, 
            row_id, 
            col_id, 
            tracer.GetRadiusPixel()
        );    
        
        // TODO: if new range is not complete => intensity=low
        if (row_id > n_row || row_id < 0 || 
            col_id > n_col || col_id < 0)
        {
            n_outrange ++;
        }
    }
    if (n_outrange == _n_cam_use)
    {
        return SMALLNUMBER;
    }

    // TODO: how to calculate ratio?
    //    1. algebraic average
    //    2. geometric average
    //    3. harmonic average 
    //    4. current one 
    std::vector<double> peak_intensity(_n_cam_use, 0); // local peak intensity 
    std::vector<std::vector<double>> otf_param_list(_n_cam_use);

    for (int m = 0; m < _n_cam_use; m ++)
    {
        otf_param_list[m] = _otf.GetOTFParam(
            _cam_id[m],
            pos_new
        );
    }

    // Updating particle intensity 
    // (ignoring the camera with highest local peak intensity)
    // Calculate intensity for each camera: 
    //  remove low intensity according to median and then remove outlier
    //  less than two cameras -> ghost, intensity = 0.
    double numerator   = 0.0;     
    double denominator = 0.0;
    std::vector<double> numerator_list(_n_cam_use, 0);
    std::vector<double> denominator_list(_n_cam_use, 0);
    for (int m = 0; m < _n_cam_use; m ++)
    {
        int row_min = search_range_list_new[m]._row_min;
        int row_max = search_range_list_new[m]._row_max;
        int col_min = search_range_list_new[m]._col_min;
        int col_max = search_range_list_new[m]._col_max;

        int row_min_old = search_range_list[m]._row_min;
        int row_max_old = search_range_list[m]._row_max;
        int col_min_old = search_range_list[m]._col_min;
        int col_max_old = search_range_list[m]._col_max;

        for (int row_id = row_min; row_id < row_max; row_id ++)
        {
            for (int col_id = col_min; col_id < col_max; col_id ++)
            {
                if ( row_id >= row_min_old && row_id < row_max_old && 
                     col_id >= col_min_old && col_id < col_max_old )
                {
                    int i = row_id - row_min_old;
                    int j = col_id - col_min_old;

                    numerator_list[m] += aug_img_list[m](i, j); 
                }
                else 
                {
                    numerator_list[m] += _res_img_list[_cam_id[m]](row_id, col_id);   
                }
                
                double value = GaussianProjection(
                    tracer.GetMatchPosInfo()[m], // how to store Match pos info?
                    otf_param_list[m],
                    col_id, row_id
                );
                denominator_list[m] += value;

                // debug
                if (value == 0)
                {
                    std::cout << "GaussianProjection: " << value << std::endl;
                    std::cout << _n_cam << "," << _n_cam_use << ",";
                    std::cout << col_id << "," << row_id << ",";
                    std::cout << tracer_center_pixel(0,0) << "," << tracer_center_pixel(1,0) << "\n"; 
                }
            }
        }
    }

    // Remove smallest proj intensity
    double median = myMATH::Median(denominator_list);
    std::vector<int> cam_id_list;
    std::vector<double> denominator_list_new;
    for (int m = 0; m < _n_cam_use; m ++)
    {
        if (denominator_list[m] >= median * min_ratio_of_median) // use a const threshold or coeff?
        {
            cam_id_list.push_back(m);
            denominator_list_new.push_back(denominator_list[m]);
        }
    }


    // Find outlier 
    std::vector<int> ignore_cam_id_list(cam_id_list.size(), false);
    myMATH::IsOutlier<double> (ignore_cam_id_list, denominator_list_new);
    int n_normal = 0;
    for (int i = 0; i < cam_id_list.size(); i ++)
    {
        if (! ignore_cam_id_list[i]) // if not outlier
        {
            int cam_id = cam_id_list[i];
            numerator += numerator_list[cam_id];
            denominator += denominator_list[cam_id];
            n_normal ++;
        }
    }
    if (n_normal < 2)
    {
        numerator = 0;
    }
              
    // Is it reasonable to calculate the ratio in this way?
    if (denominator == 0)
    {
        std::cerr << "Warning: Shake<T>::CalObjectIntensity: " 
                  << "denominator is 0!"
                  << std::endl;
        
        denominator = INTSMALLNUMBER;
    }

    double intensity_new = intensity * std::sqrt(
        std::fabs(numerator/denominator)
    );
    
    return intensity_new;
}


template<class T>
double Shake<T>::ShakingTracer (TracerInfo& tracer, double delta, double intensity)
{
    // std::cout << "qhh debug -1" << std::endl;

    std::vector<PixelRange> search_range_list(_n_cam_use);
    std::vector<Matrix<double>> aug_img_list;
    // for (int i = 0; i < _n_cam_use; i ++)
    // {
    //     aug_img_list.push_back(Matrix<double>(tracer.GetRadiusPixel()*2+2, tracer.GetRadiusPixel()*2+2));
    // }

    // std::cout << "qhh debug 0" << std::endl;

    Matrix<double> tracer_center_pixel(3,1,0);
    for (int m = 0; m < _n_cam_use; m ++)
    {
        // Get the actual pixel range
        tracer_center_pixel = tracer.GetMatchPosInfo()[m];
        
        PixelRange search_range = FindSearchRegion(
            m, 
            tracer_center_pixel(1,0), 
            tracer_center_pixel(0,0), 
            tracer.GetRadiusPixel()
        );
        
        // std::cout << "qhh debug 0.1" << std::endl;

        #pragma region AugmentedImage 
        // TODO: if 2d search range is not full => directly set intensity as 0
        // Creating a particle reproj image (I_p) matrix in the pixel range
        if (search_range.GetNumOfRow() <= 0 || 
            search_range.GetNumOfCol() <= 0)
        {
            std::cout << "qhh debug 0.2" 
                      << "n_row = " << search_range.GetNumOfRow() << ", "
                      << "n_col = " << search_range.GetNumOfCol()
                      << std::endl;
            std::cout << "tracer_pixel = " << std::endl;
            tracer_center_pixel.Print();
            tracer.GetCenterPos().Print();
            continue;
        }

        Matrix<double> aug_img(
            search_range.GetNumOfRow(),
            search_range.GetNumOfCol()
        );

        // std::cout << "qhh debug 0.3" << std::endl;

        std::vector<double> otf_para = _otf.GetOTFParam(_cam_id[m], tracer.GetCenterPos());

        int i = 0;
        for (int row_id = search_range._row_min; row_id < search_range._row_max; row_id++)
        {
            int j = 0;
            for (int col_id = search_range._col_min; col_id < search_range._col_max; col_id++)
            {
                double value = GaussianProjection(
                    tracer_center_pixel,
                    otf_para,
                    col_id, row_id
                );

                // Creating a particle augmented residual image: (I_(res+p))
                // TODO: try 1.5 search range
                aug_img(i, j) = std::max( 0.0,
                    std::min(
                        double(_max_intensity[_cam_id[m]]),
                        _res_img_list[_cam_id[m]](row_id, col_id) + value
                    )
                );
                j++;
            }
            i++;
        }
        #pragma endregion create augmented image

        // std::cout << "qhh debug 0.4" << std::endl;

        search_range_list[m] = search_range;
        aug_img_list.push_back(aug_img);

    }

    // std::cout << "qhh debug 1" << std::endl;

    // Updating the particle position and intensity
    Matrix<double> pos_new (
        UpdatePos3D (
            tracer.GetCenterPos(), 
            search_range_list, 
            aug_img_list, 
            delta,
            tracer.GetRadiusPixel()
        )
    );
    tracer.SetCenterPos(pos_new);

    // std::cout << "qhh debug 2" << std::endl;

    // Update 2d match position 
    tracer.ClearMatchPosInfo();
    tracer.SetMatchPosInfo(_cam_use);

    // std::cout << "qhh debug 3" << std::endl;

    // Updating the tracer intensity
    // sum up the intensity of all cam 
    // to delete ghost tracer 
    double object_intensity = CalObjectIntensity(
        tracer,
        search_range_list,
        aug_img_list, 
        intensity,
        _int_low_ratio
    );

    // std::cout << "qhh debug 4" << std::endl;

    return object_intensity;
}


template<class T>
void Shake<T>::TracerResImg ()
{
    Matrix<double> pt_world(3,1,0);
    Matrix<double> pt_img(3,1,0);
    int n_object = _object_info.size();
    for (int i = 0; i < n_object; i ++)
    {
        pt_world = _object_info[i].GetCenterPos();

        for (int m = 0; m < _n_cam_use; m ++)
        {
            std::vector<double> otf_para = _otf.GetOTFParam(_cam_id[m], pt_world);
            pt_img = _object_info[i].GetMatchPosInfo(m);
            PixelRange search_range = FindSearchRegion(
                m, 
                pt_img(1,0),
                pt_img(0,0),
                _object_info[i].GetRadiusPixel()
            );

            int id = _cam_id[m];
            int row_min = search_range._row_min;
            int row_max = search_range._row_max;
            int col_min = search_range._col_min;
            int col_max = search_range._col_max;
            for (int row_id = row_min; row_id < row_max; row_id ++)
            {
                for (int col_id = col_min; col_id < col_max; col_id ++)
                {
                    // _res_img_list[id](row_id, col_id) -= 
                    //     GaussianProjection(
                    //         pt_img,
                    //         otf_para,
                    //         col_id, row_id
                    // );
                    _res_img_list[id](row_id, col_id) = _orig_img_list[id](row_id, col_id) 
                        - GaussianProjection(
                            pt_img,
                            otf_para,
                            col_id, row_id
                    );

                    if (!std::isfinite(_res_img_list[id](row_id, col_id)))
                    {
                        std::cerr << "Shake<T>::TracerResImg pos2: "
                                << "cam_id = " << id << ", "
                                << "(row_id, col_id) = " 
                                << "(" << row_id << ","
                                << col_id << "), "
                                << "_res_img_list[cam_id](row_id, col_id) = "
                                << _res_img_list[id](row_id, col_id) << "\n"
                                << "otf_para = " << otf_para[0] << "," 
                                                 << otf_para[1] << ","
                                                 << otf_para[2] << ","
                                                 << otf_para[3] << "."
                                                 << "\n"
                                << std::endl;
                        std::cerr << "pt_img = " << std::endl;
                        pt_img.Print();
                        std::cerr << "pt_world = " << std::endl;
                        pt_world.Print();
                        throw;
                    }
                }
            }
        }
    }    
}


template<class T>
void Shake<T>::AbsResImg ()
{
    int id;
    for (int m = 0; m < _n_cam_use; m ++)
    {
        id = _cam_id[m];
        for (int row_id = 0; row_id < _res_img_list[id].GetDimX(); row_id ++)
        {
            for (int col_id = 0; col_id < _res_img_list[id].GetDimY(); col_id ++)
            {
                if (!std::isfinite(_res_img_list[id](row_id, col_id)))
                {
                    std::cerr << "Shake<T>::AbsResImg: "
                              << "cam_id = " << id << ", "
                              << "(row_id, col_id) = " 
                              << "(" << row_id << ","
                              << col_id << "), "
                              << "_res_img_list[cam_id](row_id, col_id) = "
                              << _res_img_list[id](row_id, col_id)
                              << std::endl;
                    throw;
                }

                _res_img_list[id](row_id, col_id) = std::max(
                    0.0, 
                    _res_img_list[id](row_id, col_id)
                );
            }
        }
    }
}


template<class T>
void Shake<T>::RemoveGhost ()
{
    int n_object = _object_info.size();

    // compute mean intensity for the smaller part 
    double mean_int = 0; // mean_int_temp = 0;

    int count = 0;
    std::vector<int> outlier(n_object, false);
    myMATH::IsOutlier(outlier, _intensity_list);
    for (int i = 0; i < n_object; i ++)
    {
        if (! outlier[i])
        {
            mean_int += _intensity_list[i];
            count ++;
        }   
    }
    mean_int /= count;

    // Remove ghost tracers based on intensity 
    for (int i = n_object-1; i >= 0; i --)
    {
        if (_intensity_list[i] < _int_low_ratio * mean_int) // use a const threshold or coeff?
        {
            _intensity_list.erase(_intensity_list.begin() + i);
            _remove_object_id.push_back(i);
            _object_info.erase(_object_info.begin() + i);
            _n_ghost ++;
        }
        else
        {
            _keep_object_id_back.push_back(i);
        }
    }
}


template<class T>
void Shake<T>::RunShake (bool TRIONLY)
{ 
    int n_object = _object_info.size();

    // For tracers
    if (typeid(T) == typeid(TracerInfo))
    {
        if (!TRIONLY)
        {
            // Initialize residue img 
            double delta;
            for (int loop = 0; loop < _n_loop; loop ++)
            {
                // std::cout << "Shake loop " << loop << " start!" << std::endl;

                // for (int cam_id = 0; cam_id < _n_cam; cam_id ++)
                // {
                //     _res_img_list[cam_id] = _orig_img_list[cam_id];
                // }
                TracerResImg();

                if (loop < 1)
                {
                    delta = _shake_width;
                }
                else if (loop < 5)
                {
                    delta = _shake_width / std::pow(2, loop - 1);   
                }
                else 
                {
                    delta = _shake_width / 20;
                }

                #pragma omp parallel
                {
                    #pragma omp for
                    for (int i = 0; i < n_object; i ++)
                    {
                        _intensity_list[i] = ShakingTracer(_object_info[i], delta, _intensity_list[i]);
                    }
                }
            }

            // Remove ghost objects 
            RemoveGhost ();
        }
        

        // Compute final res img 
        TracerResImg();
        
        // Calculate absolute intensity value
        AbsResImg();
    }
    

}

#endif