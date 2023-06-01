#include <typeinfo>

#include <vector>
// #include <deque>
#include <typeinfo>
#include <omp.h>
#include <ctime>

#include "myMATH.h"
#include "Matrix_bool.h"
#include "Matrix.h"
#include "myIO.h"
#include "Camera.h"
#include "ObjectInfo.h"
#include "ImageIO.h"
#include "ObjectFinder.h"
#include "StereoMatch.h"
#include "OTF.h"
#include "Shake.h"
#include "IPR.h"


int main()
{
    std::cout << "Simple test!" << std::endl;

    int n_cam = 4;
    std::vector<Camera> cam_list;
    std::string path_main = "./Data/";
    for (int i = 0; i < n_cam; i ++)
    {
        Camera cam(
            path_main 
            + "cam" + std::to_string(i+1) 
            + "Params.txt"
        );
        cam_list.push_back(cam);
    }
    
    std::vector<ImageIO> img_list;
    for (int i = 0; i < n_cam; i ++)
    {
        ImageIO img;
        img.LoadImgPath(
            path_main, 
            "cam" + std::to_string(i+1) + "ImageNames.txt"
        );
        img_list.push_back(img);
    }

    // Load Img
    Matrix<double> intensity(1024, 1024);
    std::vector<Matrix<double>> orig_img_list;
    std::vector<int> max_intensity_list(n_cam, 0);
    std::vector<int> min_intensity_list(n_cam, 0);
    for (int i = 0; i < n_cam; i ++)
    {
        intensity = img_list[i].LoadImg(0);
        orig_img_list.push_back(intensity);

        max_intensity_list[i] = 255;
        min_intensity_list[i] = 10;
    }

    // OTF 
    AxisLimit limit;
    limit._x_min = -20; limit._x_max = 20;
    limit._y_min = -20; limit._y_max = 20;
    limit._z_min = -20; limit._z_max = 20;
    OTF otf(4, 3, limit);

    // IPR
    IPR<TracerInfo> ipr(orig_img_list, max_intensity_list, min_intensity_list, cam_list, otf);
    std::vector<TracerInfo> object_info;
    ipr.RunIPR(object_info);

    std::cout << "Test ended!" << std::endl;
    return 0;
}
