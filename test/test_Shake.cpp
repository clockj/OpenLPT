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

#include <time.h>
#include <random>

// test stereomatch from image with deleting ghost
bool test_function_1 ()
{
    std::cout << "test_function_1" << std::endl;

    CamList cam_list;
    for (int i = 0; i < 4; i ++)
    {
        Camera cam("../test/inputs/test_Shake/cam" + std::to_string(i+1) + ".txt");
        cam_list.cam_list.push_back(cam);
        cam_list.intensity_max.push_back(255); // default
        cam_list.useid_list.push_back(i);
    }

    // load image path
    std::vector<ImageIO> imgio_list;
    for (int i = 0; i < 4; i ++)
    {
        ImageIO imgio;
        imgio.loadImgPath("../test/inputs/test_Shake/", "cam" + std::to_string(i+1) + "ImageNames" + ".txt");
        imgio_list.push_back(imgio);
    }

    // load image
    std::vector<Image> img_list;
    for (int i = 0; i < 4; i ++)
    {
        img_list.push_back(imgio_list[i].loadImg(0));
    }

    // find 2d tracer
    std::vector<std::vector<Tracer2D>> tr2d_list_all;
    std::vector<double> properties = {255, 30, 2};
    ObjectFinder2D objfinder;
    for (int i = 0; i < 4; i ++)
    {
        std::vector<Tracer2D> tr2d_list;
        objfinder.findObject2D(tr2d_list, img_list[i], properties);
        tr2d_list_all.push_back(tr2d_list);
        std::cout << "tr2d_list_all[" << i << "].size() = " << tr2d_list_all[i].size() << std::endl;
    }

    // stereo match
    SMParam param;
    param.tor_2d = 1.;
    param.tor_3d = 2.4e-2;
    param.n_thread = 6;
    param.check_id = 3;
    param.check_radius = 3; // try 4
    param.is_delete_ghost = true;
    param.is_update_inner_var = false;
    StereoMatch stereo_match(param, cam_list);

    clock_t start, end;
    start = clock();
    std::vector<Tracer3D> tr3d_list;
    stereo_match.match(tr3d_list, tr2d_list_all);
    end = clock();
    std::cout << "match time = " << double(end-start)/CLOCKS_PER_SEC << " [s]" << std::endl;

    // save tracer info
    stereo_match.saveObjInfo("../test/results/test_Shake/tr3d_img_deleteGhost.csv", tr3d_list);

    // Shake
    std::vector<Tracer3D> tr3d_list_shake(tr3d_list);

    // add noise
    std::default_random_engine generator(1234);
    std::normal_distribution<double> dist(0, 0.001);
    for (int i = 0; i < tr3d_list_shake.size(); i ++)
    {
        tr3d_list_shake[i]._pt_center[0] += dist(generator);
        tr3d_list_shake[i]._pt_center[1] += dist(generator);
        tr3d_list_shake[i]._pt_center[2] += dist(generator);
    }

    AxisLimit boundary;
    boundary.x_min = -20;
    boundary.x_max = 20;
    boundary.y_min = -20;
    boundary.y_max = 20;
    boundary.z_min = -20;
    boundary.z_max = 20;
    OTF otf;
    otf.loadParam(4, 2, 2, 2, boundary);

    // Shake s (cam_list, 0.01, 0.1, 4, 6); // 0.25 vox, 1 vox = 0.04 mm
    Shake s (cam_list, 0.01, 50, 10, 6); // 0.25 vox, 1 vox = 0.04 mm
    // Shake s (cam_list, 2*param.tor_3d, 3, 4, 6); // 0.25 vox, 1 vox = 0.04 mm

    // Shake s (cam_list, 0.01, 1000, 1, 6); // gradient descent

    start = clock();
    s.runShake(tr3d_list_shake, otf, img_list, false); // 0.25 vox, 1 vox = 0.04 mm
    end = clock();
    std::cout << "shake time = " << double(end-start)/CLOCKS_PER_SEC << " [s]" << std::endl;
    
    // save tracer info after shaking
    stereo_match.saveObjInfo("../test/results/test_Shake/tr3d_img_shake.csv", tr3d_list_shake);

    // load solution: pt3d_list
    Matrix<double> pt3d_list_sol("../test/solutions/test_Shake/pt3d_list_img.csv");

    // check match correctness
    int n_tr3d_real = pt3d_list_sol.getDimRow();
    int n_tr3d_find = tr3d_list.size();
    std::vector<int> is_mismatch(n_tr3d_real, 0);
    // double tor = param.tor_3d; // [mm]
    double tor = 1e-3; // [mm]
    
    // before shake
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
    std::cout << "Before shake" << std::endl;
    std::cout << "n_tr3d_real = " << n_tr3d_real << std::endl;
    std::cout << "n_tr3d_find = " << n_tr3d_find << std::endl;
    std::cout << "n_mismatch = " << n_mismatch << std::endl;
    std::cout << "notfind_rate = " << notfind_rate << "%" << std::endl;
    std::cout << "correct_rate = " << correct_rate << "%" << std::endl;

    // after shake
    std::fill(is_mismatch.begin(), is_mismatch.end(), 0);
    n_tr3d_find = tr3d_list_shake.size();
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
            error_list[j] = myMATH::dist(pt3d, tr3d_list_shake[j]._pt_center);
        }
                
        std::vector<double>::iterator minIt = std::min_element(error_list.begin(), error_list.end());

        if (*minIt > tor)
        {
            is_mismatch[i] = 1;
        }
    }
    
    n_mismatch = 0;
    for (int i = 0; i < n_tr3d_real; i ++)
    {
        n_mismatch += is_mismatch[i];
    }

    notfind_rate = double(n_mismatch) / n_tr3d_real * 100; // [%]   
    correct_rate = double(n_tr3d_real-n_mismatch) / n_tr3d_find * 100; // [%]
    std::cout << "After shake" << std::endl;
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
    fs::create_directory("../test/results/test_Shake/");

    IS_TRUE(test_function_1());

    return 0;
}