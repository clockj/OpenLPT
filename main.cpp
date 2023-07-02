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
#include "Track.h"
#include "STB.h"


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
    std::vector<std::vector<TracerInfo>> pt_list_all;
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
        pt_list_all.push_back(object_info);
    }

    // Predicted field
    std::cout << "PredField start!" << std::endl;

    // // for debug
    // Matrix<double> mtx_pre("3Dpoints_pf_1.csv");
    // std::vector<Matrix<double>> pt_list_pre(mtx_pre.GetDimX(), Matrix<double>(3,1,0));
    // for (int i = 0; i < mtx_pre.GetDimX(); i ++)
    // {
    //     pt_list_pre[i](0,0) = mtx_pre(i,0);
    //     pt_list_pre[i](1,0) = mtx_pre(i,1);
    //     pt_list_pre[i](2,0) = mtx_pre(i,2);
    // }

    // Matrix<double> mtx_cur("3Dpoints_pf_2.csv");
    // std::vector<Matrix<double>> pt_list_cur(mtx_cur.GetDimX(), Matrix<double>(3,1,0));
    // for (int i = 0; i < mtx_cur.GetDimX(); i ++)
    // {
    //     pt_list_cur[i](0,0) = mtx_cur(i,0);
    //     pt_list_cur[i](1,0) = mtx_cur(i,1);
    //     pt_list_cur[i](2,0) = mtx_cur(i,2);
    // }

    // pt_list_all.push_back(pt_list_pre);
    // pt_list_all.push_back(pt_list_cur);


    std::vector<int> n_xyz = {50,50,50};
    double r = 25 * 40/1000; // 25 * (max-min)/1000
    PredField<TracerInfo> pf (limit, n_xyz, pt_list_all[0], pt_list_all[1], r);

    Matrix<double> disp_pred(3,1);
    disp_pred = pf.PtInterp(pt_list_all[1][0].GetCenterPos());
    // disp_pred.Print();
    Matrix<double> disp_field(3,50*50*50);
    disp_field = pf.GetField();
    disp_field.WriteMatrix("pf_disp_field.csv");

    STB<TracerInfo> stb("./Data/stbConfig.txt");

    std::cout << "Finish OpenLPT!" << std::endl;
    return 0;
}
