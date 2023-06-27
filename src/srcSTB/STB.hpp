#ifndef STB_HPP
#define STB_HPP

#include "STB.h"


template<class T>
STB<T>::STB (std::string folder)
{
    // Read param
    std::ifstream stb_config(folder + "stbConfig.txt", std::ios::in);
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

    
    // Load _cam_list & _imgio_list
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


    // 


    // Load particle radius
    parsed >> _r_object;
    

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


    // 

}




#endif