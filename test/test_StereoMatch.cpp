#include "test.h"
#include "Matrix.h"
#include "Camera.h"
#include "ImageIO.h"
#include "ObjectInfo.h"
#include "ObjectFinder.h"
#include "StereoMatch.h"

#include <time.h>

// test tracer stereomatch without deleting ghost
bool test_function_1 ()
{    
    std::cout << "test_function_1" << std::endl;

    CamList cam_list;
    for (int i = 0; i < 4; i ++)
    {
        Camera cam("../test/inputs/test_StereoMatch/cam" + std::to_string(i+1) + ".txt");
        cam_list.cam_list.push_back(cam);
        cam_list.intensity_max.push_back(255); // default
        cam_list.useid_list.push_back(i);
    }

    // load 2d tracer
    std::vector<std::vector<Tracer2D>> tr2d_list_all;
    Tracer2D tr2d;
    for (int i = 0; i < 4; i ++)
    {
        std::vector<Tracer2D> tr2d_list;
        
        std::ifstream file("../test/inputs/test_StereoMatch/pt2d_list_cam" + std::to_string(i+1) + ".csv");
        std::string line;
        while (std::getline(file, line))
        {
            std::stringstream ss(line);
            std::string token;
            std::vector<std::string> tokens;
            while (std::getline(ss, token, ','))
            {
                tokens.push_back(token);
            }
            tr2d._pt_center[0] = std::stod(tokens[0]);
            tr2d._pt_center[1] = std::stod(tokens[1]);
            tr2d_list.push_back(tr2d);
        }
        file.close();

        tr2d_list_all.push_back(tr2d_list);
    }

    // stereo match
    SMParam param;
    param.tor_2d = 1;
    param.tor_3d = 5e-3;
    param.n_thread = 6;
    param.check_id = 3;
    param.check_radius = 1;
    param.is_delete_ghost = false;
    StereoMatch stereo_match(param, cam_list);

    clock_t start, end;
    start = clock();
    std::vector<Tracer3D> tr3d_list;
    stereo_match.match(tr3d_list, tr2d_list_all);
    end = clock();
    std::cout << "time = " << double(end-start)/CLOCKS_PER_SEC << " [s]" << std::endl;

    // save tracer info
    stereo_match.saveObjInfo("../test/results/test_StereoMatch/tr3d.csv", tr3d_list);
    stereo_match.saveObjIDMatchList("../test/results/test_StereoMatch/objID_match_list.csv");

    // load solution: pt3d_list
    Matrix<double> pt3d_list_sol("../test/solutions/test_StereoMatch/pt3d_list.csv");

    // check match correctness
    int n_tracer = pt3d_list_sol.getDimRow();
    int n_match = stereo_match._objID_match_list.size();

    int id;
    bool is_correct;
    std::vector<int> correctID_list;
    for (int i = 0; i < n_match; i ++)
    {
        id = stereo_match._objID_match_list[i][0];
        is_correct = true;
        for (int j = 0; j < stereo_match._n_cam_use; j ++)
        {
            if (stereo_match._objID_match_list[i][j] != id)
            {
                is_correct = false;
                break;
            }
        }
        if (is_correct)
        {
            correctID_list.push_back(i);
        }
    }

    int n_correct = correctID_list.size();
    int n_ghost = n_match - n_correct;
    double correct_rate = (double)n_correct / n_tracer * 100; // [%]
    double ghost_rate = (double)n_ghost / n_tracer * 100; // [%]

    std::cout << "test_function_1:" << std::endl;
    std::cout << "n_tracer = " << n_tracer << std::endl;
    std::cout << "n_match = " << n_match << std::endl;
    std::cout << "n_correct = " << n_correct << std::endl;
    std::cout << "n_ghost = " << n_ghost << std::endl;
    std::cout << "correct% = " << correct_rate << "%" << std::endl;
    std::cout << "ghost% = " << ghost_rate << "%" << std::endl;

    double correct_rate_thres = 100;
    double ghost_rate_thres = 5;
    if (correct_rate < correct_rate_thres || ghost_rate > ghost_rate_thres)
    {
        std::cout << "test_function_1 error at line: " << __LINE__ << std::endl;
        std::cout << "correct_rate < " << correct_rate_thres << "%" 
                  << "or ghost_rate > " << ghost_rate_thres << "%" << std::endl;
        return false;
    }

    // check correctness of 3D position with correct match
    double tor = 1e-5; // [mm]
    double error;
    int tr_id, sol_id;
    Pt3D pt3d_sol;
    for (int i = 0; i < n_correct; i ++)
    {
        tr_id = correctID_list[i];
        sol_id = stereo_match._objID_match_list[tr_id][0];

        pt3d_sol[0] = pt3d_list_sol(sol_id, 0);
        pt3d_sol[1] = pt3d_list_sol(sol_id, 1);
        pt3d_sol[2] = pt3d_list_sol(sol_id, 2);

        error = myMATH::dist(pt3d_sol, tr3d_list[tr_id]._pt_center);

        if (error > tor)
        {
            std::cout << "test_function_1 error at line: " << __LINE__ << std::endl;
            std::cout << "error = " << error << " > " << tor << std::endl;
            std::cout << "tr_id = " << tr_id << std::endl;
            std::cout << "sol_id = " << sol_id << std::endl;
            std::cout << "pt3d_sol = ";
            pt3d_sol.print();
            std::cout << "tr3d_list[tr_id]._pt_center = ";
            tr3d_list[tr_id]._pt_center.print();

            return false;
        }
    }

    std::cout << "test_function_1 passed\n" << std::endl;

    return true;
}

// test tracer stereomatch with deleting ghost
bool test_function_2 ()
{
    std::cout << "test_function_2" << std::endl;

    CamList cam_list;
    for (int i = 0; i < 4; i ++)
    {
        Camera cam("../test/inputs/test_StereoMatch/cam" + std::to_string(i+1) + ".txt");
        cam_list.cam_list.push_back(cam);
        cam_list.intensity_max.push_back(255); // default
        cam_list.useid_list.push_back(i);
    }

    // load 2d tracer
    std::vector<std::vector<Tracer2D>> tr2d_list_all;
    Tracer2D tr2d;
    for (int i = 0; i < 4; i ++)
    {
        std::vector<Tracer2D> tr2d_list;
        
        std::ifstream file("../test/inputs/test_StereoMatch/pt2d_list_cam" + std::to_string(i+1) + ".csv");
        std::string line;
        while (std::getline(file, line))
        {
            std::stringstream ss(line);
            std::string token;
            std::vector<std::string> tokens;
            while (std::getline(ss, token, ','))
            {
                tokens.push_back(token);
            }
            tr2d._pt_center[0] = std::stod(tokens[0]);
            tr2d._pt_center[1] = std::stod(tokens[1]);
            tr2d_list.push_back(tr2d);
        }
        file.close();

        tr2d_list_all.push_back(tr2d_list);
    }

    // stereo match
    SMParam param;
    param.tor_2d = 1;
    param.tor_3d = 5e-3;
    param.n_thread = 6;
    param.check_id = 3;
    param.check_radius = 1;
    param.is_delete_ghost = true;
    param.is_update_inner_var = true;
    StereoMatch stereo_match(param, cam_list);

    clock_t start, end;
    start = clock();
    std::vector<Tracer3D> tr3d_list;
    stereo_match.match(tr3d_list, tr2d_list_all);
    end = clock();
    std::cout << "time = " << double(end-start)/CLOCKS_PER_SEC << " [s]" << std::endl;

    // save tracer info
    stereo_match.saveObjInfo("../test/results/test_StereoMatch/tr3d_deleteGhost.csv", tr3d_list);
    stereo_match.saveObjIDMatchList("../test/results/test_StereoMatch/objID_match_list_deleteGhost.csv");

    // load solution: pt3d_list
    Matrix<double> pt3d_list_sol("../test/solutions/test_StereoMatch/pt3d_list.csv");

    // check match correctness
    int n_tracer = pt3d_list_sol.getDimRow();
    int n_match = stereo_match._objID_match_list.size();

    int id;
    bool is_correct;
    std::vector<int> correctID_list;
    for (int i = 0; i < n_match; i ++)
    {
        id = stereo_match._objID_match_list[i][0];
        is_correct = true;
        for (int j = 0; j < stereo_match._n_cam_use; j ++)
        {
            if (stereo_match._objID_match_list[i][j] != id)
            {
                is_correct = false;
                break;
            }
        }
        if (is_correct)
        {
            correctID_list.push_back(i);
        }
    }

    int n_correct = correctID_list.size();
    int n_ghost = n_match - n_correct;
    double correct_rate = (double)n_correct / n_tracer * 100; // [%]
    double ghost_rate = (double)n_ghost / n_tracer * 100; // [%]

    std::cout << "test_function_1:" << std::endl;
    std::cout << "n_tracer = " << n_tracer << std::endl;
    std::cout << "n_match = " << n_match << std::endl;
    std::cout << "n_correct = " << n_correct << std::endl;
    std::cout << "n_ghost = " << n_ghost << std::endl;
    std::cout << "correct% = " << correct_rate << "%" << std::endl;
    std::cout << "ghost% = " << ghost_rate << "%" << std::endl;

    double correct_rate_thres = 99;
    double ghost_rate_thres = 5;
    if (correct_rate < correct_rate_thres || ghost_rate > ghost_rate_thres)
    {
        std::cout << "test_function_2 error at line: " << __LINE__ << std::endl;
        std::cout << "correct_rate < " << correct_rate_thres << "%" 
                  << "or ghost_rate > " << ghost_rate_thres << "%" << std::endl;
        return false;
    }

    // check correctness of 3D position with correct match
    double tor = 1e-5; // [mm]
    double error;
    int tr_id, sol_id;
    Pt3D pt3d_sol;
    for (int i = 0; i < n_correct; i ++)
    {
        tr_id = correctID_list[i];
        sol_id = stereo_match._objID_match_list[tr_id][0];

        pt3d_sol[0] = pt3d_list_sol(sol_id, 0);
        pt3d_sol[1] = pt3d_list_sol(sol_id, 1);
        pt3d_sol[2] = pt3d_list_sol(sol_id, 2);

        error = myMATH::dist(pt3d_sol, tr3d_list[tr_id]._pt_center);

        if (error > tor)
        {
            std::cout << "test_function_2 error at line: " << __LINE__ << std::endl;
            std::cout << "error = " << error << " > " << tor << std::endl;
            std::cout << "tr_id = " << tr_id << std::endl;
            std::cout << "sol_id = " << sol_id << std::endl;
            std::cout << "pt3d_sol = ";
            pt3d_sol.print();
            std::cout << "tr3d_list[tr_id]._pt_center = ";
            tr3d_list[tr_id]._pt_center.print();

            return false;
        }
    }

    std::cout << "test_function_2 passed\n" << std::endl;

    return true;
}

// test stereomatch from image with deleting ghost
bool test_function_3 ()
{
    std::cout << "test_function_3" << std::endl;

    CamList cam_list;
    for (int i = 0; i < 4; i ++)
    {
        Camera cam("../test/inputs/test_StereoMatch/cam" + std::to_string(i+1) + ".txt");
        cam_list.cam_list.push_back(cam);
        cam_list.intensity_max.push_back(255); // default
        cam_list.useid_list.push_back(i);
    }

    // load image path
    std::vector<ImageIO> imgio_list;
    for (int i = 0; i < 4; i ++)
    {
        ImageIO imgio;
        imgio.loadImgPath("../test/inputs/test_StereoMatch/", "cam" + std::to_string(i+1) + "ImageNames" + ".txt");
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
    param.check_radius = 3;
    param.is_delete_ghost = true;
    param.is_update_inner_var = false;
    StereoMatch stereo_match(param, cam_list);

    clock_t start, end;
    start = clock();
    std::vector<Tracer3D> tr3d_list;
    stereo_match.match(tr3d_list, tr2d_list_all);
    end = clock();
    std::cout << "time = " << double(end-start)/CLOCKS_PER_SEC << " [s]" << std::endl;

    // save tracer info
    stereo_match.saveObjInfo("../test/results/test_StereoMatch/tr3d_img_deleteGhost.csv", tr3d_list);

    // load solution: pt3d_list
    Matrix<double> pt3d_list_sol("../test/solutions/test_StereoMatch/pt3d_list_img.csv");

    // check match correctness
    int n_tr3d_real = pt3d_list_sol.getDimRow();
    int n_tr3d_find = tr3d_list.size();
    std::vector<int> is_mismatch(n_tr3d_real, 0);
    double tor = param.tor_3d; // [mm]
    
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

    double mismatch_rate = double(n_mismatch) / n_tr3d_real * 100; // [%]
    std::cout << "n_tr3d_real = " << n_tr3d_real << std::endl;
    std::cout << "n_tr3d_find = " << n_tr3d_find << std::endl;
    std::cout << "n_mismatch = " << n_mismatch << std::endl;
    std::cout << "mismatch_rate = " << mismatch_rate << "%" << std::endl;

    std::cout << "test_function_3 passed\n" << std::endl;

    return true;
}


int main ()
{
    fs::create_directories("../test/results/test_StereoMatch/");

    IS_TRUE(test_function_1());
    IS_TRUE(test_function_2());
    IS_TRUE(test_function_3());

    return 0;
}