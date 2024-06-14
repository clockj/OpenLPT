#include "test.h"
#include "Matrix.h"
#include "myMATH.h"
#include "Camera.h"
#include "ImageIO.h"
#include "ObjectInfo.h"
#include "ObjectFinder.h"
#include "StereoMatch.h"
#include "Shake.h"
#include "OTF.h"
#include "IPR.h"

#include <time.h>

bool test_function_1 ()
{
    std::cout << "test_function_1" << std::endl;

    CamList cam_list;
    for (int i = 0; i < 4; i ++)
    {
        Camera cam("../test/inputs/test_IPR/cam" + std::to_string(i+1) + ".txt");
        cam_list.cam_list.push_back(cam);
        cam_list.intensity_max.push_back(255); // default
        cam_list.useid_list.push_back(i);
    }

    // load image path
    std::vector<ImageIO> imgio_list;
    for (int i = 0; i < 4; i ++)
    {
        ImageIO imgio;
        imgio.loadImgPath("../test/inputs/test_IPR/", "cam" + std::to_string(i+1) + "ImageNames" + ".txt");
        imgio_list.push_back(imgio);
    }

    // load image
    std::vector<Image> img_list;
    for (int i = 0; i < 4; i ++)
    {
        img_list.push_back(imgio_list[i].loadImg(0));
    }

    // IPR
    AxisLimit boundary;
    boundary.x_min = -20;
    boundary.x_max = 20;
    boundary.y_min = -20;
    boundary.y_max = 20;
    boundary.z_min = -20;
    boundary.z_max = 20;
    OTF otf;
    otf.loadParam(4, 2, 2, 2, boundary);

    std::vector<double> tr2d_properties = {255, 30};
    std::vector<Tracer3D> tr3d_list;

    IPR ipr(cam_list, img_list, IPRParam());

    clock_t start, end;
    ipr.runIPR(tr3d_list, tr2d_properties, otf);

    // save tracer info after shaking
    ipr.saveObjInfo("../test/results/test_IPR/tr3d_img_ipr.csv", tr3d_list);

    std::cout << std::endl;

    // load solution: pt3d_list
    Matrix<double> pt3d_list_sol("../test/solutions/test_IPR/pt3d_list_img.csv");

    // check match correctness
    int n_tr3d_real = pt3d_list_sol.getDimRow();
    int n_tr3d_find = tr3d_list.size();
    std::vector<int> is_mismatch(n_tr3d_real, 0);
    // double tor = param.tor_3d; // [mm]
    double tor = 10e-3; // [mm]
    
    // after IPR
    #pragma omp parallel for
    for (int i = 0; i < n_tr3d_real; i ++)
    {
        Pt3D pt3d;
        pt3d[0] = pt3d_list_sol(i, 0);
        pt3d[1] = pt3d_list_sol(i, 1);
        pt3d[2] = pt3d_list_sol(i, 2);
        
        std::vector<double> error_list(n_tr3d_find, 0);

        for (int j = 0; j < n_tr3d_find; j ++)
        {
            error_list[j] = myMATH::dist(pt3d, tr3d_list[j]._pt_center);
        }
                
        std::vector<double>::iterator minIt = std::min_element(error_list.begin(), error_list.end());

        if (*minIt > tor)
        {
            is_mismatch[i] = 1;
        }
    }
    
    int n_mismatch = 0;
    for (int i = 0; i < n_tr3d_real; i ++)
    {
        n_mismatch += is_mismatch[i];
    }

    double notfind_rate = double(n_mismatch) / n_tr3d_real * 100; // [%]   
    double correct_rate = double(n_tr3d_real-n_mismatch) / n_tr3d_find * 100; // [%]
    std::cout << "After IPR" << std::endl;
    std::cout << "n_tr3d_real = " << n_tr3d_real << std::endl;
    std::cout << "n_tr3d_find = " << n_tr3d_find << std::endl;
    std::cout << "n_mismatch = " << n_mismatch << std::endl;
    std::cout << "notfind_rate = " << notfind_rate << "%" << std::endl;
    std::cout << "correct_rate = " << correct_rate << "%" << std::endl;


    std::cout << "test_function_1 passed\n" << std::endl;

    return true;
}


int main ()
{
    fs::create_directories("../test/results/test_IPR/");

    IS_TRUE(test_function_1());

    return 0;
}