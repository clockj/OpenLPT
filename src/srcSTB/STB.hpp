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
    for (int i = 0; i < n_cam; i ++)
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
    parsed >> _ipr_flag; // 1 for using ipr in initialphase (0 for use .mat files)
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
                  << "is not included in StereoMatch!"
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
        parsed >>  _imgint_max_list[i];
    }
    for (int i = 0; i < _n_cam; i ++)
    {
        parsed >>  _imgint_min_list[i];
    }
}




#endif