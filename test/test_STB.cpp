#include "test.h"

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <variant>

#include "Matrix.h"
#include "Camera.h"
#include "ObjectInfo.h"
#include "STB.h"

bool test_function_1 ()
{
    std::cout << "test_function_1" << std::endl;

    // Load config
    int frame_start, frame_end, n_thread;
    CamList cam_list;
    int n_cam_all = 4;
    AxisLimit axis_limit;
    double vx_to_mm;

    std::ifstream stb_config("../test/inputs/test_STB/config.txt", std::ios::in);
    std::vector<std::string> lines;
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
        }
        else if (commentpos == 0)
        {
            continue;
        }

        lines.push_back(line);
    }
    stb_config.close();

    std::stringstream parsed;

    // Load frame range
    int line_id = 0;
    parsed.str(lines[line_id]);
    std::getline(parsed, line, ',');
    frame_start = std::stoi(line);
    std::getline(parsed, line, ',');
    frame_end = std::stoi(line);
    parsed.clear();
    
    // Load frame rate
    line_id ++;
    int fps = std::stoi(lines[line_id]);

    // Load thread number
    line_id ++;
    n_thread = std::stoi(lines[line_id]);

    // Load cam files
    line_id ++;
    n_cam_all = std::stoi(lines[line_id]);
    for (int i = 0; i < n_cam_all; i++)
    {
        line_id ++;
        parsed.str(lines[line_id]);
        
        std::getline(parsed, line, ',');
        cam_list.cam_list.push_back(Camera(line));

        std::getline(parsed, line, ',');
        cam_list.intensity_max.push_back(std::stoi(line));

        cam_list.useid_list.push_back(i);

        parsed.clear();
    }

    // Load image io
    std::vector<ImageIO> imgio_list;
    for (int i = 0; i < n_cam_all; i++)
    {
        line_id ++;
        imgio_list.push_back(ImageIO());
        imgio_list[i].loadImgPath("", lines[line_id]);
    }

    // Load axis limit
    line_id ++;
    parsed.str(lines[line_id]);
    std::getline(parsed, line, ',');
    axis_limit.x_min = std::stod(line);
    std::getline(parsed, line, ',');
    axis_limit.x_max = std::stod(line);
    std::getline(parsed, line, ',');
    axis_limit.y_min = std::stod(line);
    std::getline(parsed, line, ',');
    axis_limit.y_max = std::stod(line);
    std::getline(parsed, line, ',');
    axis_limit.z_min = std::stod(line);
    std::getline(parsed, line, ',');
    axis_limit.z_max = std::stod(line);
    parsed.clear();

    // Load vx_to_mm
    line_id ++;
    vx_to_mm = std::stod(lines[line_id]);

    // Load output folder path
    line_id ++;
    std::string output_folder = lines[line_id];

    // Load STB
    line_id ++;    
    parsed.str(lines[line_id]);
    std::vector<std::variant<STB<Tracer3D>>> stb_list;
    int n_tr_class = 0;
    while (std::getline(parsed, line, ','))
    {
        line_id ++;
        if (line == "Tracer")
        {
            stb_list.push_back(STB<Tracer3D>(frame_start, frame_end, 1, vx_to_mm, n_thread, output_folder+"Tracer_"+std::to_string(n_tr_class)+'/', cam_list, axis_limit, lines[line_id]));

            // Calibrate OTF // 
            std::cout << "Start Calibrating OTF!" << std::endl;
            int n_otf_calib = 1;
            int n_obj2d_max = 1000;
            // int r_otf_calib = 2; // [px]
            std::vector<Image> img_list(n_otf_calib);

            std::visit(
                [&](auto& stb) 
                { 
                    double r_otf_calib = stb.getObjParam()[2];
                    for (int i = 0; i < n_cam_all; i ++)
                    {
                        for (int j = 0; j < n_otf_calib; j ++)
                        {
                            img_list[j] = imgio_list[i].loadImg(frame_start + j);
                        }
                        stb.calibrateOTF(i, n_obj2d_max, r_otf_calib, img_list); 
                    }
                }, 
                stb_list[n_tr_class]
            );

            std::cout << "Finish Calibrating OTF!\n" << std::endl;

            n_tr_class ++;
        }
        else
        {
            std::cout << "Error: Unknown object type: " << line << std::endl;
            return false;
        }
    }
    parsed.clear();

    // Start processing
    std::vector<Image> img_list(n_cam_all);
    for (int frame_id = frame_start; frame_id < frame_end+1; frame_id ++)
    {
        // load image
        for (int i = 0; i < n_cam_all; i ++)
        {
            img_list[i] = imgio_list[i].loadImg(frame_id);
        }

        for (int i = 0; i < stb_list.size(); i ++)
        {
            std::visit(
                [&](auto& stb) 
                { 
                    stb.processFrame(frame_id, img_list); 
                }, 
                stb_list[i]
            );
        }
    }

    std::cout << std::endl;

    return true;
}

bool test_function_2 ()
{
    std::cout << "test_function_2" << std::endl;

    // Load config
    int frame_start, frame_end, n_thread;
    CamList cam_list;
    int n_cam_all;
    AxisLimit axis_limit;
    double vx_to_mm;

    std::ifstream stb_config("../test/inputs/test_STB/config.txt", std::ios::in);
    std::vector<std::string> lines;
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
        }
        else if (commentpos == 0)
        {
            continue;
        }

        lines.push_back(line);
    }
    stb_config.close();

    std::stringstream parsed;

    // Load frame range
    int line_id = 0;
    parsed.str(lines[line_id]);
    std::getline(parsed, line, ',');
    frame_start = std::stoi(line);
    std::getline(parsed, line, ',');
    frame_end = std::stoi(line);
    parsed.clear();
    
    // Load frame rate
    line_id ++;
    int fps = std::stoi(lines[line_id]);

    // Load thread number
    line_id ++;
    n_thread = std::stoi(lines[line_id]);

    // Load cam files
    line_id ++;
    n_cam_all = std::stoi(lines[line_id]);
    for (int i = 0; i < n_cam_all; i++)
    {
        line_id ++;
        parsed.str(lines[line_id]);
        
        std::getline(parsed, line, ',');
        cam_list.cam_list.push_back(Camera(line));

        std::getline(parsed, line, ',');
        cam_list.intensity_max.push_back(std::stoi(line));

        cam_list.useid_list.push_back(i);

        parsed.clear();
    }

    // Load image io
    std::vector<ImageIO> imgio_list;
    for (int i = 0; i < n_cam_all; i++)
    {
        line_id ++;
        imgio_list.push_back(ImageIO());
        imgio_list[i].loadImgPath("", lines[line_id]);
    }

    // Load axis limit
    line_id ++;
    parsed.str(lines[line_id]);
    std::getline(parsed, line, ',');
    axis_limit.x_min = std::stod(line);
    std::getline(parsed, line, ',');
    axis_limit.x_max = std::stod(line);
    std::getline(parsed, line, ',');
    axis_limit.y_min = std::stod(line);
    std::getline(parsed, line, ',');
    axis_limit.y_max = std::stod(line);
    std::getline(parsed, line, ',');
    axis_limit.z_min = std::stod(line);
    std::getline(parsed, line, ',');
    axis_limit.z_max = std::stod(line);
    parsed.clear();

    // Load vx_to_mm
    line_id ++;
    vx_to_mm = std::stod(lines[line_id]);

    // Load output folder path
    line_id ++;
    std::string output_folder = lines[line_id];

    // Calibrate OTF
    int n_otf_calib = 1;
    int n_obj2d_max = 1000;
    int r_otf_calib = 2; // [px]
    OTF otf(n_cam_all, 2, 2, 2, axis_limit);
    std::vector<Image> img_list(n_otf_calib);
    
    for (int i = 0; i < n_cam_all; i ++)
    {
        std::vector<std::vector<Tracer2D>> tr2d_list_all;
        std::cout << "Camera " << i << std::endl;
        std::cout << "\tNumber of found tracers in each image: ";
        for (int j = 0; j < n_otf_calib; j ++)
        {
            img_list[j] = imgio_list[i].loadImg(frame_start + j);

            ObjectFinder2D objfinder;
            std::vector<Tracer2D> tr2d_list;
            objfinder.findObject2D(tr2d_list, img_list[j], {double(cam_list.intensity_max[i]), 30., 2.});

            std::cout << tr2d_list.size();

            // if tr2d_list is too large, randomly select some tracers
            int seed = 123;
            if (tr2d_list.size() > n_obj2d_max)
            {
                std::shuffle(tr2d_list.begin(), tr2d_list.end(), std::default_random_engine(seed));
                tr2d_list.erase(tr2d_list.begin()+n_obj2d_max, tr2d_list.end());

                std::cout << "(" << n_obj2d_max << ")";
            }
            else if (tr2d_list.size() == 0)
            {
                std::cerr << "Quit OTF calibration: No tracer found in camera " << i << std::endl;
                throw error_size;
            }

            std::cout << ",";

            tr2d_list_all.push_back(tr2d_list);
        }

        std::cout << std::endl;

        calibOTFParam(otf._param, i, r_otf_calib, tr2d_list_all, img_list);
    }
    
    std::cout << "a = " << std::endl;
    otf._param.a.print();
    std::cout << "b = " << std::endl;
    otf._param.b.print();
    std::cout << "c = " << std::endl;
    otf._param.c.print();
    std::cout << "alpha = " << std::endl;
    otf._param.alpha.print();

    // Calculate residue image
    int frame_id = 0;

    #pragma omp parallel for
    for (int i = 0; i < n_cam_all; i ++)
    {
        double a, b, c, alpha, xc, yc, dx, dy, intensity;
        int row_min, row_max, col_min, col_max, n_row, n_col;
        Image img_res, img_orig;
        ObjectFinder2D objfinder;

        img_orig = imgio_list[i].loadImg(frame_id);

        std::vector<Tracer2D> tr2d_list;
        objfinder.findObject2D(tr2d_list, img_orig, {double(cam_list.intensity_max[i]), 30., 2.});

        // calculate residue image
        a = otf._param.a(i,0);
        b = otf._param.b(i,0);
        c = otf._param.c(i,0);
        n_row = img_orig.getDimRow();
        n_col = img_orig.getDimCol();

        img_res = img_orig;
        for (int j = 0; j < tr2d_list.size(); j ++)
        {
            xc = tr2d_list[j]._pt_center[0];
            yc = tr2d_list[j]._pt_center[1];

            row_min = std::max(0, int(yc - tr2d_list[j]._r_px));
            row_max = std::min(n_row, int(yc + tr2d_list[j]._r_px + 1));
            col_min = std::max(0, int(xc - tr2d_list[j]._r_px));
            col_max = std::min(n_col, int(xc + tr2d_list[j]._r_px + 1));

            for (int row = row_min; row < row_max; row ++)
            {
                for (int col = col_min; col < col_max; col ++)
                {
                    dx = col - xc;
                    dy = row - yc;

                    intensity = a * std::exp(- b * dx*dx - c * dy*dy);
                    intensity = img_orig(row, col) - intensity;

                    img_res(row, col) = intensity < img_res(row,col) ? intensity : img_res(row,col);
                }
            }
        }

        img_res.write("../test/results/test_STB/"+std::to_string(frame_id)+"_img_res_"+std::to_string(i)+".csv");
    }


    return true;
}

int main()
{
    fs::create_directories("../test/results/test_STB/");

    IS_TRUE(test_function_1());
    // IS_TRUE(test_function_2());

    return 0;
}