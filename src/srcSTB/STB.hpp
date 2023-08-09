#ifndef STB_HPP
#define STB_HPP

#include "STB.h"


template<class T>
STB<T>::STB (std::string file) :_ngrid_xyz(3,0)
{
    // Read param
    std::ifstream stb_config(file, std::ios::in);
    std::stringstream parsed;
    std::string line;
    while (std::getline(stb_config, line))
    {
        size_t commentpos = line.find('#');
        if (commentpos > 0)
        {
            if (commentpos < std::string::npos)
            {
                line.erase(commentpos);
            }
            parsed << line << '\t';
        }
    }
    stb_config.close();


    //        //
    // Camera //
    //        //
    parsed >> _n_cam;
    parsed >> _img_folder;
    for (int i = 0; i < _n_cam; i ++)
    {
        Camera cam(
            _img_folder 
            + "cam" + std::to_string(i+1) 
            + "Params.txt"
        );
        _cam_list.push_back(cam);
    }
    _n_pix_w = _cam_list[0].GetNpixw();
    _n_pix_h = _cam_list[1].GetNpixh();
    for (int i = 0; i < _n_cam; i ++)
    {
        ImageIO img;
        img.LoadImgPath(
            _img_folder, 
            "cam" + std::to_string(i+1) + "ImageNames.txt"
        );
        _imgio_list.push_back(img);
    }


    //          //
    // Tracking //
    //          //
    parsed >> _first; _first -= 1;
    parsed >> _last;  _last -= 1;

    // View area
    parsed >> _xyz_limit._x_min;
    parsed >> _xyz_limit._x_max;
    parsed >> _xyz_limit._y_min;
    parsed >> _xyz_limit._y_max;
    parsed >> _xyz_limit._z_min;
    parsed >> _xyz_limit._z_max;
    OTF otf(4, 3, _xyz_limit);

    parsed >> _vox_to_mm;

    // Initial phase
    parsed >> _ipr_flag; // 1 for using ipr in initialphase (0 for use .csv files)
    parsed >> _r_search_init; _r_search_init *= _vox_to_mm; // mm

    // Convergence phase
    parsed >> _shake_shift; _shake_shift *= _vox_to_mm;
    parsed >> _space_avg; _space_avg *= _vox_to_mm;
    parsed >> _shift_betframe_max; _shift_betframe_max *= _vox_to_mm;
    parsed >> _shift_abs_max; _shift_abs_max *= _vox_to_mm;
    parsed >> _shift_rel_max;
    parsed >> _fpt;
    parsed >> _int_lower;


    //                  //
    // Predictive Field //
    //                  //
    int n_grid;
    parsed >> n_grid; _ngrid_xyz[0] = n_grid;
    parsed >> n_grid; _ngrid_xyz[1] = n_grid;
    parsed >> n_grid; _ngrid_xyz[2] = n_grid;


    //     //
    // IPR //
    //     //
    // Load IPR option
    int option;
    parsed >> option;
    if (option == 1)
    {
        _TRIONLY = true;
    }
    else if (option == 2)
    {
        _IPRONLY = true;
    }

    // Load IPR iteration
    parsed >> _ipr_loop_outer;
    parsed >> _ipr_loop_inner;
    parsed >> _ipr_int_lower;
    parsed >> _ipr_tol_2d; _ipr_tol_2d *= _vox_to_mm;
    parsed >> _ipr_tol_3d; _ipr_tol_3d *= _vox_to_mm;

    parsed >> _ipr_is_reduced;
    parsed >> _ipr_loop_outer_reduced;
    parsed >> _ipr_tol_2d_reduced; _ipr_tol_2d_reduced *= _vox_to_mm;
    parsed >> _ipr_tol_3d_reduced; _ipr_tol_3d_reduced *= _vox_to_mm;
    
    //             //
    // Object Info //
    //             //
    if (typeid(T) == typeid(TracerInfo))
    {   
        SetTracerParam (parsed);
    }
    else 
    {
        std::cerr << "StereoMatch: " 
                  << "class " << typeid(T).name()
                  << " is not included in StereoMatch!"
                  << std::endl;
        throw;
    }
}


template<class T>
void STB<T>::SetTracerParam (std::stringstream& parsed)
{
    // Load particle radius
    int r_tracer;
    parsed >> r_tracer;
    _obj_info_sample.SetRadiusPixel(r_tracer);

    _imgint_max_list = std::vector<int>(_n_cam, 0);
    _imgint_min_list = std::vector<int>(_n_cam, 0);
    for (int i = 0; i < _n_cam; i ++)
    {
        parsed >> _imgint_max_list[i];
    }
    for (int i = 0; i < _n_cam; i ++)
    {
        parsed >> _imgint_min_list[i];
    }
}


template<class T>
void STB<T>::InitialPhase ()
{
    int endframe = _last + 1;
    if (!_TRIONLY && !_IPRONLY)
    {
        endframe = (_last >= _first + 4) ? _first+4 : _last+1;
    }

    Matrix<double> intensity(_n_pix_h, _n_pix_w);
    for (int i = 0; i < _n_cam; i ++)
    {
        _img_org_list.push_back(intensity);
        // _img_res_list.push_back(intensity);
    }

    for (int frame = _first; frame < endframe; frame ++)
    {
        std::cout << "IPR on frame " << frame << " of " << _last << std::endl;

        if (_ipr_flag)
        {
            // Load img
            for (int i = 0; i < _n_cam; i ++)
            {
                _img_org_list[i] = _imgio_list[i].LoadImg(frame);
            }

            // IPR
            IPR<T> ipr(_img_org_list, _imgint_max_list, _imgint_min_list, _cam_list, _otf);
            ipr.SetTRIONLY(_TRIONLY);
            ipr.SetIPRTimes(_ipr_loop_outer);
            ipr.SetShakeTimes(_ipr_loop_inner);
            ipr.SetTol2D(_ipr_tol_2d);
            ipr.SetTol3D(_ipr_tol_3d);

            std::vector<T> obj_list;
            ipr.RunIPR(obj_list, _ipr_is_reduced);
            _ipr_matched.push_back(obj_list);
        }
    }

    if (!_TRIONLY && !_IPRONLY)
    {
        for (int frame = _first; frame < endframe-1; frame ++)
        {   
            int currframe = frame;
            int nextframe = frame + 1;

            std::cout << "STB initial phase tracking b/w frames: "
                      << currframe << " & " << nextframe;
            
            PredField<T> pf(_xyz_limit, _ngrid_xyz, _ipr_matched[currframe-_first], _ipr_matched[nextframe-_first], _r_search_pred);

            // link particles b/w two frames
            clock_t t_start, t_end;
            t_start = clock();

            

            t_end = clock();
            std::cout << "; link particles time: " 
                      << (double) (t_end - t_start)/CLOCKS_PER_SEC
                      << std::endl;
        }
    }
    
}


// template<class T>




#endif