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
// #include "OTF.h"
// #include "Shake.h"
// #include "IPR.h"

// For test 
// std::vector<Matrix<double>> Load3DPoints (std::string file_name, int n_pt)
// {
//     std::ifstream infile(file_name, std::ios::in);
    
//     std::vector<Matrix<double>> pt_list_3d;
//     for (int i = 0; i < n_pt; i ++)
//     {
//         Matrix<double> pt(3,1);
//         infile >> pt[0] >> pt[1] >> pt[2];
//         pt_list_3d.push_back(pt);
//     }

//     return pt_list_3d;
// }


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

    Matrix<int> intensity(1024, 1024);
    ObjectFinder<TracerInfo> tf;
    int depth = 255;
    int threshold = 10;

    // Find tracer list @ t=0
    std::vector<std::vector<TracerInfo>> tracer_list_pixel;
    for (int i = 0; i < n_cam; i ++)
    {
        intensity = img_list[i].LoadImg(0);
        std::vector<TracerInfo> tracer_found 
            = tf.FindObject(intensity, depth, threshold);
        std::cout << "Number of found tracer: " << tracer_found.size() << std::endl;
        tracer_list_pixel.push_back(tracer_found);
    }
    // myIO::WriteTracerPos("Result/cam_1_tracer_2D.csv", tracer_list_pixel[0]);
    // myIO::WriteTracerPos("Result/cam_2_tracer_2D.csv", tracer_list_pixel[1]);

    // Stereomatch
    clock_t t_start, t_end;

    t_start = clock();
    StereoMatch<TracerInfo> tracer_match(
        cam_list,
        2e-2, //2e-2,
        6e-2  //1e-3
    );

    // tracer_match.SetNumThread(4);
    // tracer_match.SetTolerance1D(3);
    // tracer_match.Match(tracer_list_pixel, 0);
    tracer_match.Match(tracer_list_pixel, 1);
    
    // std::vector<std::vector<int>> tracer_match_list; 
    // tracer_match.GetObjectIDMatchList(tracer_match_list);

    std::vector<TracerInfo> tracer_list;
    tracer_match.GetObjectInfoMatchList(tracer_list);

    std::vector<std::vector<int>> match_list;
    tracer_match.GetObjectIDMatchList(match_list);

    std::vector<double> error_list;
    tracer_match.GetTriErrorList(error_list);
    
    t_end = clock();

    std::cout << "Matching time: " 
              << (double) (t_end - t_start)/CLOCKS_PER_SEC
              << std::endl;

    std::cout << "Test ended!" << std::endl;
    return 0;
}
