#ifndef IPR_HPP
#define IPR_HPP

#include "IPR.h"

template<class T>
IPR<T>::IPR(std::vector<Matrix<double>>& orig_img_list, std::vector<int>& max_intensity_list, std::vector<int>& min_intensity_list, std::vector<Camera>& cam_list, OTF& otf) 
    : _orig_img_list(orig_img_list), _max_intensity_list(max_intensity_list), _min_intensity_list(min_intensity_list), _cam_list(cam_list), _n_cam(cam_list.size()), _otf(otf), _res_img_list(orig_img_list) 
{}

template<class T>
void IPR<T>::RunIPR(std::vector<T>& object_info, bool is_reduced, int n_reduced) 
{
    std::vector<Matrix<double>> intensity_list(_orig_img_list);

    for (int loop_id = 0; loop_id < _n_loop_ipr; loop_id ++)
    {
        // Step 1: identify 2D position from image
        std::vector<std::vector<T>> tracer_list_pixel;
        for (int i = 0; i < _n_cam; i ++)
        {
            ObjectFinder<T> tf;
            std::vector<T> tracer_found 
                = tf.FindObject(_res_img_list[i], _max_intensity_list[i], _min_intensity_list[i]);
            std::cout << "\t Number of found tracer: " << tracer_found.size() << std::endl;

            tracer_list_pixel.push_back(tracer_found);
        }


        // Step 2: stereomatching
        clock_t t_start, t_end;
        t_start = clock();

        StereoMatch<T> tracer_match(
            _cam_list,
            _tol_2d,
            _tol_3d 
        );
        tracer_match.Match(tracer_list_pixel, 1);
        
        std::vector<T> tracer_info_match_list;
        tracer_match.GetObjectInfoMatchList(
            tracer_info_match_list
        );

        t_end = clock();
        std::cout << "\t Matching time: " 
                  << (double) (t_end - t_start)/CLOCKS_PER_SEC
                  << std::endl;


        // Step 3: Shaking
        if (tracer_info_match_list.size() == 0)
        {
            break;
        }
        // std::cout << "Start shake!" << std::endl;
        Shake<T> s(
            _res_img_list,
            tracer_info_match_list, // also used as output
            _shake_width, // unit: mm
            _cam_list,
            _otf
        );
        s.SetShakeTimes (_n_loop_shake);
        s.RunShake();
        // std::cout << "Finish shake!" << std::endl;

        
        s.GetResImg(_res_img_list);
        

        _object_info.insert(_object_info.end(), tracer_info_match_list.begin(), tracer_info_match_list.end());

        std::cout << "\t IPR step " << loop_id << ": find " << _object_info.size() << " particles. " << std::endl;
    }


    // Reduced cam
    if (is_reduced)
    {
        int n_rest = _n_cam - n_reduced;
        if (_n_cam - n_reduced < 2 || n_reduced < 0)
        {
            std::cerr << "IPR RunIPR() problem! "
                      << "n_cam-n_reduced = " << n_rest
                      << "out of range."
                      << std::endl;
            throw;
        }

        std::cout << "\t Start reduced camera loop!" << std::endl;

        // Generate all possible cam_id
        std::deque<std::vector<int>> cam_id_all;
        std::vector<int> cam_id;
        CreateCamID (cam_id_all, cam_id, 0, n_rest);

        for (int i = 0; i < cam_id_all.size(); i ++)
        {
            ReducedCamLoop(cam_id_all[i], _n_loop_ipr_reduced);
        }
    }


    std::cout << "IPR Finish! " << "Find " << _object_info.size() << " particles. " << std::endl;
    object_info = _object_info;
}

// pending 07/04
template<class T>
void IPR<T>::ReducedCamLoop(std::vector<int> cam_id, int n_loop)
{
    int n_cam = cam_id.size();
    int id = 0;
    std::vector<Camera> cam_list;
    for (int i = 0; i < n_cam; i ++)
    {
        id = cam_id[i];
        cam_list.push_back(_cam_list[id]);
    }

    for (int loop = 0; loop < n_loop; loop ++)
    {
        // Step 1: identify 2D position from image
        std::vector<std::vector<T>> tracer_list_pixel;
        for (int i = 0; i < n_cam; i ++)
        {
            id = cam_id[i];
            ObjectFinder<T> tf;
            std::vector<T> tracer_found 
                = tf.FindObject(_res_img_list[id], _max_intensity_list[id], _min_intensity_list[id]);
            std::cout << "\t Number of found tracer: " << tracer_found.size() << std::endl;

            tracer_list_pixel.push_back(tracer_found);
        }

        // Step 2: stereomatching
        clock_t t_start, t_end;
        t_start = clock();

        StereoMatch<T> tracer_match(
            cam_list,
            _tol_2d,
            _tol_3d 
        );
        tracer_match.Match(tracer_list_pixel, 1);
        
        std::vector<T> tracer_info_match_list;
        tracer_match.GetObjectInfoMatchList(
            tracer_info_match_list
        );

        t_end = clock();
        std::cout << "\t Matching time: " 
                    << (double) (t_end - t_start)/CLOCKS_PER_SEC
                    << std::endl;

        // Step 3: Shaking
        if (tracer_info_match_list.size() == 0)
        {
            break;
        }
        Shake<T> s(
            _res_img_list,
            tracer_info_match_list, // also used as output
            _shake_width, // unit: mm
            _cam_list,
            _otf,
            cam_id
        );
        s.SetShakeTimes (_n_loop_shake);
        s.RunShake();
  
        s.GetResImg(_res_img_list);
        
        _object_info.insert(_object_info.end(), tracer_info_match_list.begin(), tracer_info_match_list.end());
    }
    
}


template<class T>
void IPR<T>::CreateCamID (std::deque<std::vector<int>>& cam_id_all, std::vector<int> cam_id, int id, int n_rest)
{
    int n_cam_id = cam_id.size();
    if (n_cam_id < n_rest && id < _n_cam-(n_rest-n_cam_id)+1)
    {
        for (int i = id; i < _n_cam-(n_rest-n_cam_id)+1; i ++)
        {
            cam_id.push_back(i);
            CreateCamID(cam_id_all, cam_id, i+1, n_rest);
            cam_id.pop_back();
        }
    }
    else
    {
        if (n_cam_id == n_rest)
        {
            cam_id_all.push_back(cam_id);
        }
        else
        {
            std::cerr << "IPR CreateCamID problem!" 
                      << "cam_id.size=" << cam_id.size()
                      << ", n_rest=" << n_rest 
                      << std::endl;
            throw;
        }
    }
}

#endif