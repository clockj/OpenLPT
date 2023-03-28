#ifndef SHAKE_HPP
#define SHAKE_HPP

#include "Shake.h"

template<class T>
PixelRange Shake<T>::FindSearchRegion (int cam_id, int center_row, int center_col, int half_width_pixel)
{
    // int center_row = tracer.tracer_center_pixel(1,0);
    // int center_col = tracer.tracer_center_pixel(0,0);

    PixelRange search_range;
    // why max(1,...) in the original code?
    // but with min(Npixw,...)
    // I choose to start from 0
    search_range._row_min = std::max(
        0,
        center_row - half_width_pixel
    );
    search_range._col_min = std::max(
        0,
        center_col - half_width_pixel
    );
    search_range._row_max = std::min(
        _cam[cam_id].GetNpixh(),
        center_row + half_width_pixel
    );
    search_range._col_max = std::min(
        _cam[cam_id].GetNpixw(),
        center_col + half_width_pixel
    );

    return search_range;
}


// (x, y) = (col_id, row_id)
template<class T>
double Shake<T>::GaussianProjection (Matrix<double>& tracer_center_pixel, std::vector<double>& otf_param, int x, int y)
{
    double xx =   (double(x) - tracer_center_pixel(0,0)) * std::cos(otf_param[3]) 
                + (double(y) - tracer_center_pixel(1,0)) * std::sin(otf_param[3]);
    double yy = - (double(x) - tracer_center_pixel(0,0)) * std::sin(otf_param[3]) 
                + (double(y) - tracer_center_pixel(1,0)) * std::cos(otf_param[3]);
    double value = otf_param[0] 
                   * std::exp( - (otf_param[1] * xx*xx + otf_param[2] * yy*yy) );
    return value;
}


template<class T>
double Shake<T>::CalPointResidue (Matrix<double> pos_new, std::vector<PixelRange>& search_range_list, std::vector<Matrix<double>>& aug_img_list)
{
    double residue = 0;
    int n_cam = _cam.size();

    for (int cam_id = 0; cam_id < n_cam; cam_id ++)
    {
        Matrix<double> pos_2d_pixel(
            _cam[cam_id].ImgMMToPixel(
                _cam[cam_id].WorldToImgMM(pos_new)
            )
        );

        std::vector<double> otf_param = _otf.GetOTFParam(cam_id, pos_new);

        // Calculating the residual for updated 3D position
        // use old range
        // how to define residue?
        // in the original code, 
        //  residue = sum (sum residue^2) for all cam 
        // should it be 
        //  residue = sum sqrt(sum residue^2) / n_cam ?
        // I still use the original definition here 
        PixelRange search_range = search_range_list[cam_id];
        int i = 0; 
        for (int row_id = search_range._row_min; row_id < search_range._row_max; row_id ++)
        {
            int j = 0; 
            for (int col_id = search_range._col_min; col_id < search_range._col_max; col_id ++)
            {
                residue += std::pow(
                    aug_img_list[cam_id](i,j)
                    - std::round(
                        GaussianProjection(
                            pos_2d_pixel,
                            otf_param,
                            col_id, row_id
                        )
                    ),
                    2
                );
                j ++;
            }
            i ++;
        }
    }

    return residue;
}


// only use the initial 3 elements
template<class T>
std::vector<double> Shake<T>::FitQuadratic (double* array, std::vector<double>& residue)
{
    std::vector<double> coeff(3,0);

    // a 
    coeff[0] = (
                (residue[1]-residue[0]) / (array[1]-array[0]) 
                - (residue[2]-residue[0]) / (array[2]-array[0])
               ) / (array[1]-array[2]);
    // b 
    coeff[1] = (residue[1]-residue[0]) / (array[1]-array[0]) 
                - coeff[0] * (array[1]+array[0]);
    // c 
    coeff[2] = residue[2] - coeff[0]*array[2]*array[2] - coeff[1]*array[2];

    return coeff;
}


template<class T>
Matrix<double> Shake<T>::UpdatePos3D (Matrix<double>& pos_old, std::vector<PixelRange>& search_range_list, std::vector<Matrix<double>>& aug_img_list, double delta)
{
    Matrix<double> pos_new(pos_old);

    double delta_list[3];
    delta_list[0] = - delta;
    delta_list[1] = 0;
    delta_list[2] = delta;

    double array_list[4] = {0};

    // shaking on x,y,z direction
    for (int i = 0; i < 3; i ++)
    {
        std::vector<double> coeff;
        std::vector<double> residue_list;

        for (int j = 0; j < 3; j ++)
        {
            array_list[j] = pos_old(i,0) + delta_list[j];
            pos_new(i,0) = array_list[j];
            residue_list.push_back(
                CalPointResidue(
                    pos_new,
                    search_range_list,
                    aug_img_list
                )
            );
        }
        coeff = FitQuadratic(array_list, residue_list);

        array_list[3] = - coeff[1] / (2 * coeff[0]);
        if (array_list[3]>array_list[0] && array_list[3]<array_list[2])
        {
            pos_new(i,0) = array_list[3];
            residue_list.push_back(
                CalPointResidue(
                    pos_new,
                    search_range_list,
                    aug_img_list
                )
            );
        }
        pos_new(i,0) = array_list[myMATH::BubbleSort(residue_list)[0]];
        // TODO: update Augemented image according to new position?
    }

    return pos_new;
}


// remove search_range_small, directly use new range 
// if use new range, if it is out of the original range, use the original image intensity
template<class T>
double Shake<T>::CalObjectIntensity (TracerInfo& tracer, std::vector<PixelRange>& search_range_list, std::vector<Matrix<double>>& aug_img_list, double intensity, double min_ratio_of_median)
{
    int n_cam = _cam.size();
    Matrix<double> pos_new = tracer.GetCenterPos();

    // Pixel range on the cameras around the updated particle centers、
    std::vector<PixelRange> search_range_list_new(n_cam);
    for (int cam_id = 0; cam_id < n_cam; cam_id ++)
    {
        Matrix<double> tracer_center_pixel 
            = _cam[cam_id].ImgMMToPixel(
                _cam[cam_id].WorldToImgMM(
                    pos_new
                )
            );

        search_range_list_new[cam_id] = FindSearchRegion(
            cam_id, 
            tracer_center_pixel(1,0), 
            tracer_center_pixel(0,0), 
            tracer.GetRadiusPixel()
        );      
    }

    // TODO: how to calculate ratio?
    //    1. algebraic average
    //    2. geometric average
    //    3. harmonic average 
    //    4. current one 
    std::vector<double> peak_intensity(n_cam, 0); // local peak intensity 
    std::vector<std::vector<double>> otf_param_list(n_cam);

    for (int cam_id = 0; cam_id < n_cam; cam_id ++)
    {
        otf_param_list[cam_id] = _otf.GetOTFParam(
            cam_id,
            pos_new
        );
    }

    //     int row_min = search_range_list_new[cam_id]._row_min;
    //     int row_max = search_range_list_new[cam_id]._row_max;
    //     int col_min = search_range_list_new[cam_id]._col_min;
    //     int col_max = search_range_list_new[cam_id]._col_max;

    //     for (int row_id = row_min; row_id < row_max; row_id ++)
    //     {
    //         for (int col_id = col_min; col_id < col_max; col_id ++)
    //         {
    //             if ( row_id >= row_min && row_id < row_max && 
    //                  col_id >= col_min && col_id < col_max )
    //             {
    //                 int i = row_id - row_min;
    //                 int j = col_id - col_min;
    //                 peak_intensity[cam_id] = std::max(
    //                     peak_intensity[cam_id],
    //                     aug_img_list[cam_id](i, j)
    //                 );
    //             }
    //             else 
    //             {
    //                 peak_intensity[cam_id] = std::max(
    //                     peak_intensity[cam_id],
    //                     _res_img_list[cam_id](row_id, col_id)
    //                 );
    //             }
    //         }
    //     }
    // }
    // int ignore_cam_id_new = myMATH::MaxIndexBubbleSort(peak_intensity);
    // std::vector<int> ignore_cam_id_list = myMATH::IsOutlier(peak_intensity);


    // Updating particle intensity 
    // (ignoring the camera with highest local peak intensity)
    // Calculate intensity for each camera: 
    //  remove low intensity according to median and then remove outlier
    //  less then two cameras -> ghost, intensity = 0.
    double numerator   = 0.0;     
    double denominator = 0.0;
    std::vector<double> numerator_list(n_cam, 0);
    std::vector<double> denominator_list(n_cam, 0);
    for (int cam_id = 0; cam_id < n_cam; cam_id ++)
    {
        int row_min = search_range_list_new[cam_id]._row_min;
        int row_max = search_range_list_new[cam_id]._row_max;
        int col_min = search_range_list_new[cam_id]._col_min;
        int col_max = search_range_list_new[cam_id]._col_max;

        for (int row_id = row_min; row_id < row_max; row_id ++)
        {
            for (int col_id = col_min; col_id < col_max; col_id ++)
            {
                if ( row_id >= row_min && row_id < row_max && 
                     col_id >= col_min && col_id < col_max )
                {
                    int i = row_id - row_min;
                    int j = col_id - col_min;

                    // numerator += aug_img_list[cam_id](i, j);
                    // denominator += std::round(
                    //     GaussianProjection(
                    //         tracer.GetMatchPosInfo()[cam_id],
                    //         otf_param_list[cam_id],
                    //         col_id, row_id
                    //     )
                    // );
                    numerator_list[cam_id] += aug_img_list[cam_id](i, j); 
                }
                else 
                {
                    numerator_list[cam_id] += _res_img_list[cam_id](row_id, col_id);   
                }

                denominator_list[cam_id] += std::round(
                    GaussianProjection(
                        tracer.GetMatchPosInfo()[cam_id],
                        otf_param_list[cam_id],
                        col_id, row_id
                    )
                );
            }
        }
    }

    // Remove smallest proj intensity
    double median = myMATH::Median(denominator_list);
    std::vector<int> cam_id_list;
    std::vector<double> denominator_list_new;
    for (int cam_id = 0; cam_id < n_cam; cam_id ++)
    {
        if (denominator_list[cam_id] > median * min_ratio_of_median)
        {
            cam_id_list.push_back(cam_id);
            denominator_list_new.push_back(denominator_list[cam_id]);
        }
    }


    // Find outlier 
    std::vector<int> ignore_cam_id_list = myMATH::IsOutlier(denominator_list_new);
    int n_normal = 0;
    for (int i = 0; i < cam_id_list.size(); i ++)
    {
        if (ignore_cam_id_list[i] == 0)
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
    int n_cam = _cam.size();

    // std::cout << "qhh debug -1" << std::endl;

    std::vector<PixelRange> search_range_list(n_cam);
    std::vector<Matrix<double>> aug_img_list;
    for (int i = 0; i < n_cam; i ++)
    {
        aug_img_list.push_back(Matrix<double>(tracer.GetRadiusPixel() * 2+2, tracer.GetRadiusPixel() * 2+2));
    }

    // std::cout << "qhh debug 0" << std::endl;

    for (int cam_id = 0; cam_id < n_cam; cam_id ++)
    {
        // Get the actual pixel range
        Matrix<double> tracer_center_pixel = tracer.GetMatchPosInfo()[cam_id];
        
        PixelRange search_range = FindSearchRegion(
            cam_id, 
            tracer_center_pixel(1,0), 
            tracer_center_pixel(0,0), 
            tracer.GetRadiusPixel() * 2
        );
        
        // std::cout << "qhh debug 0.1" << std::endl;

        #pragma region AugmentedImage 
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
            continue;
        }

        Matrix<double> aug_img(
            search_range.GetNumOfRow(),
            search_range.GetNumOfCol()
        );

        // std::cout << "qhh debug 0.3" << std::endl;

        std::vector<double> otf_para = _otf.GetOTFParam(cam_id, tracer.GetCenterPos());

        int i = 0;
        for (int row_id = search_range._row_min; row_id < search_range._row_max; row_id++)
        {
            int j = 0;
            for (int col_id = search_range._col_min; col_id < search_range._col_max; col_id++)
            {
                double value = std::round(
                    GaussianProjection(
                        tracer_center_pixel,
                        otf_para,
                        col_id, row_id
                    )
                );


                // Creating a particle augmented residual image: (I_(res+p))
                // aug_img(i, j) = std::max(
                //     0.0,
                //     std::min(
                //         double(_max_intensity[cam_id]),
                //         _res_img_list[cam_id](row_id, col_id) + value
                //     )
                // );
                aug_img(i, j) = std::min(
                    double(_max_intensity[cam_id]),
                    _res_img_list[cam_id](row_id, col_id) + value
                );
                j++;
            }
            i++;
        }
        #pragma endregion create augmented image

        // std::cout << "qhh debug 0.4" << std::endl;

        search_range_list[cam_id] = search_range;

        // std::cout << "qhh debug 0.5" << std::endl;

        aug_img_list[cam_id] = aug_img;

    }

    // std::cout << "qhh debug 1" << std::endl;

    // Updating the particle position and intensity
    Matrix<double> pos_old = tracer.GetCenterPos();
    Matrix<double> pos_new = UpdatePos3D (
        pos_old, 
        search_range_list, 
        aug_img_list, 
        delta
    );
    tracer.SetCenterPos(pos_new);

    // std::cout << "qhh debug 2" << std::endl;

    // Update 2d match position 
    tracer.ClearMatchPosInfo();
    tracer.SetMatchPosInfo(_cam);

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
    int n_cam = _cam.size();

    int n_object = _object_info.size();
    for (int i = 0; i < n_object; i ++)
    {
        Matrix<double> pt_world = _object_info[i].GetCenterPos();

        for (int cam_id = 0; cam_id < n_cam; cam_id ++)
        {
            std::vector<double> otf_para = _otf.GetOTFParam(cam_id, pt_world);
            Matrix<double> pt_img = _object_info[i].GetMatchPosInfo(cam_id);
            PixelRange search_range = FindSearchRegion(
                cam_id, 
                pt_img(1,0),
                pt_img(0,0),
                _object_info[i].GetRadiusPixel() * 2
            );

            int row_min = search_range._row_min;
            int row_max = search_range._row_max;
            int col_min = search_range._col_min;
            int col_max = search_range._col_max;
            for (int row_id = row_min; row_id < row_max; row_id ++)
            {
                for (int col_id = col_min; col_id < col_max; col_id ++)
                {
                    _res_img_list[cam_id](row_id, col_id) -= 
                        GaussianProjection(
                            pt_img,
                            otf_para,
                            col_id, row_id
                    );

                    if (!std::isfinite(_res_img_list[cam_id](row_id, col_id)))
                    {
                        std::cerr << "Shake<T>::TracerResImg pos2: "
                                << "cam_id = " << cam_id << ", "
                                << "(row_id, col_id) = " 
                                << "(" << row_id << ","
                                << col_id << "), "
                                << "_res_img_list[cam_id](row_id, col_id) = "
                                << _res_img_list[cam_id](row_id, col_id) << "\n"
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
    int n_cam = _cam.size();

    for (int cam_id = 0; cam_id < n_cam; cam_id ++)
    {
        for (int row_id = 0; row_id < _res_img_list[cam_id].GetDimX(); row_id ++)
        {
            for (int col_id = 0; col_id < _res_img_list[cam_id].GetDimY(); col_id ++)
            {
                if (!std::isfinite(_res_img_list[cam_id](row_id, col_id)))
                {
                    std::cerr << "Shake<T>::AbsResImg: "
                              << "cam_id = " << cam_id << ", "
                              << "(row_id, col_id) = " 
                              << "(" << row_id << ","
                              << col_id << "), "
                              << "_res_img_list[cam_id](row_id, col_id) = "
                              << _res_img_list[cam_id](row_id, col_id)
                              << std::endl;
                    throw;
                }

                _res_img_list[cam_id](row_id, col_id) = std::max(
                    0.0, 
                    _res_img_list[cam_id](row_id, col_id)
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
    
    std::vector<int> outlier = myMATH::IsOutlier(_intensity_list);
    for (int i = 0; i < n_object; i ++)
    {
        if (!outlier[i])
        {
            mean_int += _intensity_list[i];
            count ++;
        }   
    }

    //mean_int = mean_int_temp;
    //mean_int_temp /= n_object;

    //for (int i = 0; i < n_object; i ++)
    //{
    //    if (_intensity_list[i] > 30 * mean_int_temp)
    //    {
    //        mean_int -= _intensity_list[i];
    //        count ++;
    //    }
    //}

    mean_int /= count;


    // Remove ghost tracers based on intensity 
    for (int i = n_object-1; i >= 0; i --)
    {
        if (_intensity_list[i] < _int_low_ratio * mean_int)
        {
            _intensity_list.erase(_intensity_list.begin() + i);
            _remove_object_id.push_back(i);
            _object_info.erase(_object_info.begin() + i);
            _n_ghost ++;
        }
    }
}


template<class T>
void Shake<T>::RunShake ()
{ 
    int n_object = _object_info.size();

    if (typeid(T) == typeid(TracerInfo))
    {
        // For tracers
        // Initialize intensity estimation 
        for (int i = 0; i < n_object; i ++)
        {
            _intensity_list.push_back(1.0);
        }
        // Initialize residue img 
        int n_cam = _cam.size();
        for (int cam_id = 0; cam_id < n_cam; cam_id ++)
        {
            _res_img_list.push_back(_orig_img_list[cam_id].TypeToDouble());
        }


        for (int loop = 0; loop < _n_loop; loop ++)
        {
            std::cout << "Shake loop " << loop << " start!" << std::endl;

            double delta;

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
            
            TracerResImg();

            // std::cout << "zsj shake 2" << std::endl;

            #pragma omp parallel
            {
                #pragma omp for
                for (int i = 0; i < n_object; i ++)
                {
                    // std::cout << "zsj shake 2.1, i = " << i << std::endl;
                    _intensity_list[i] = ShakingTracer(_object_info[i], delta, _intensity_list[i]);
                    // std::cout << "zsj shake 2.1, i = " << i << std::endl;
                }
            }
            
            // std::cout << "zsj shake 3" << std::endl;
        }

        // Remove ghost objects 
        RemoveGhost ();

        // std::cout << "zsj shake 4" << std::endl;

        // Compute final res img 
        TracerResImg();

        // std::cout << "zsj shake 5" << std::endl;

        AbsResImg();
    }
    

}

#endif