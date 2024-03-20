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


    return true;
}

int main()
{
    IS_TRUE(test_function_1());

    return 0;
}