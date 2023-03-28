#ifndef IPR_HPP
#define IPR_HPP

#include "IPR.h"

template<class T>
IPR<T>::IPR(std::vector<Matrix<int>>& orig_img_list, std::vector<int>& max_intensity_list, std::vector<int>& min_intensity_list, std::vector<Camera>& cam_list, OTF& otf) 
    : _orig_img_list(orig_img_list), _max_intensity_list(max_intensity_list), _min_intensity_list(min_intensity_list), _cam_list(cam_list), _otf(otf) 
{}

template<class T>
std::vector<T> IPR<T>::RunIPR() 
{
    std::vector<Matrix<int>> intensity_list(_orig_img_list);

    int n_cam = _cam_list.size();
    for (int loop_id = 0; loop_id < _n_loop_ipr; loop_id ++)
    {
        // Step 1: identify 2D position from image
        std::vector<std::vector<T>> tracer_list_pixel;
        for (int i = 0; i < n_cam; i ++)
        {
            ObjectFinder<T> tf;
            std::vector<T> tracer_found 
                = tf.FindObject(intensity_list[i], _max_intensity_list[i], _min_intensity_list[i]);
            std::cout << "Number of found tracer: " << tracer_found.size() << std::endl;

            tracer_list_pixel.push_back(tracer_found);
        }


        // Step 2: stereomatching
        clock_t t_start, t_end;
        t_start = clock();

        StereoMatch<T> tracer_match(
            _cam_list,
            _tol_2d,//3e-2, 2e-2,
            _tol_3d //8e-2, 1e-3
        );
        tracer_match.Match(tracer_list_pixel, 1);
        
        std::vector<T> tracer_info_match_list;
        tracer_match.GetObjectInfoMatchList(
            tracer_info_match_list
        );

        t_end = clock();
        std::cout << "Matching time: " 
                  << (double) (t_end - t_start)/CLOCKS_PER_SEC
                  << std::endl;


        // Step 3: Shaking
        if (tracer_info_match_list.size() == 0)
        {
            break;
        }
        std::cout << "Start shake!" << std::endl;
        Shake<T> s(
            intensity_list,
            tracer_info_match_list, // also used as output
            _shake_width, // unit: mm
            _cam_list,
            _otf
        );
        s.SetShakeTimes (_n_loop_shake);
        s.RunShake();
        std::cout << "Finish shake!" << std::endl;

        
        intensity_list = s.GetResImg();
        

        _object_info.insert(_object_info.end(), tracer_info_match_list.begin(), tracer_info_match_list.end());

        std::cout << "IPR step " << loop_id << ": find " << _object_info.size() << " particles. " << std::endl;
    }

    std::cout << "IPR Finish! " << "Find " << _object_info.size() << " particles. " << std::endl;
    return _object_info;
}

#endif