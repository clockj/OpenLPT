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
// #include "Shake.h"
// #include "IPR.h"


int main()
{
    std::cout << "Simple test!" << std::endl;

    // int n_cam = 4;
    // std::vector<Camera> cam_list;
    // std::string path_main = "./Data/";
    // for (int i = 0; i < n_cam; i ++)
    // {
    //     Camera cam(
    //         path_main 
    //         + "cam" + std::to_string(i+1) 
    //         + "Params.txt"
    //     );
    //     cam_list.push_back(cam);
    // }
    
    // std::vector<ImageIO> img_list;
    // for (int i = 0; i < n_cam; i ++)
    // {
    //     ImageIO img;
    //     img.LoadImgPath(
    //         path_main, 
    //         "cam" + std::to_string(i+1) + "ImageNames.txt"
    //     );
    //     img_list.push_back(img);
    // }

    // Matrix<int> intensity(1024, 1024);
    // ObjectFinder<TracerInfo> tf;
    // int depth = 255;
    // int threshold = 10;

    // // Find tracer list @ t=0
    // std::vector<std::vector<TracerInfo>> tracer_list_pixel;
    // for (int i = 0; i < n_cam; i ++)
    // {
    //     intensity = img_list[i].LoadImg(0);
    //     std::vector<TracerInfo> tracer_found 
    //         = tf.FindObject(intensity, depth, threshold);
    //     std::cout << "Number of found tracer: " << tracer_found.size() << std::endl;
    //     tracer_list_pixel.push_back(tracer_found);
    // }

    // // Stereomatch
    // clock_t t_start, t_end;

    // t_start = clock();
    // StereoMatch<TracerInfo> tracer_match(
    //     cam_list,
    //     1e-2, 
    //     2.4e-2  
    // ); // original code: divide each direction into 1000 voxels
    // // tracer_match.SetNumThread(4);
    // // tracer_match.SetTolerance1D(3);
    // // tracer_match.Match(tracer_list_pixel, 0);
    // tracer_match.Match(tracer_list_pixel, 1);
    // tracer_match.WriteObjectInfoMatchList(std::string("Result/tracer_list.csv"));
    
    // t_end = clock();

    // std::cout << "Matching time: " 
    //           << (double) (t_end - t_start)/CLOCKS_PER_SEC
    //           << std::endl;

    // OTF 
    AxisLimit limit;
    limit._x_min = -20; limit._x_max = 20;
    limit._y_min = -20; limit._y_max = 20;
    limit._z_min = -20; limit._z_max = 20;
    // OTF otf(4, 3, limit);
    OTF otf(4, 3, limit, "Result/param");

    Matrix<double> pt_world(3,1,0,10,10);
    std::vector<double> param = otf.GetOTFParam(0, pt_world);
    myIO::WriteMatrix<double> ("Result/param.csv", param);

    std::cout << "Test ended!" << std::endl;
    return 0;
}
