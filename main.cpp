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
#include "PredField.h"


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

    
    // OTF 
    AxisLimit limit;
    limit._x_min = -20; limit._x_max = 20;
    limit._y_min = -20; limit._y_max = 20;
    limit._z_min = -20; limit._z_max = 20;
    OTF otf(4, 3, limit);

    // Get found pt list for initial frames
    int n_frame = 2;
    Matrix<double> intensity(1024, 1024);
    std::vector<int> max_intensity_list(n_cam, 255);
    std::vector<int> min_intensity_list(n_cam, 10);
    std::vector<std::vector<Matrix<double>>> pt_list_all;
    for (int frame_id = 0; frame_id < n_frame; frame_id ++)
    {
        // Load Img
        std::vector<Matrix<double>> orig_img_list; 
        for (int i = 0; i < n_cam; i ++)
        {
            intensity = img_list[i].LoadImg(frame_id);
            orig_img_list.push_back(intensity);
        }

        // IPR
        IPR<TracerInfo> ipr(orig_img_list, max_intensity_list, min_intensity_list, cam_list, otf);
        ipr.SetShakeTimes(4);
        ipr.SetIPRTimes(4);
        std::vector<TracerInfo> object_info;
        ipr.RunIPR(object_info);
        pt_list_all.push_back(ipr.GetPtList());
    }

    // Predicted field
    std::cout << "PredField start!" << std::endl;

    std::vector<int> n_xyz = {50,50,50};
    double r = 25 * 40/1000; // 25 * (max-min)/1000
    PredField pf (limit, n_xyz, pt_list_all[0], pt_list_all[1], r);

    Matrix<double> disp_pred(3,1);
    disp_pred = pf.PtInterp(pt_list_all[1][0]);
    disp_pred.Print();

    return 0;
}
