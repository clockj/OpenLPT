#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <variant>

#include "STBCommons.h"
#include "ImageIO.h"
#include "Matrix.h"
#include "Camera.h"
#include "ObjectInfo.h"
#include "STB.h"

void run (std::string file)
{
    std::cout << "**************" << std::endl;
    std::cout << "OpenLPT start!" << std::endl;
    std::cout << "**************" << std::endl;
    std::cout << std::endl;

    // Load config
    int frame_start, frame_end, fps, n_thread;
    CamList cam_list;
    int n_cam_all;
    AxisLimit axis_limit;
    double vx_to_mm;

    std::cout << "Load config file: " << file << std::endl;
    std::ifstream stb_config(file, std::ios::in);
    std::vector<std::string> lines;
    std::string line;
    while (std::getline(stb_config, line))
    {
        // Trim leading and trailing whitespace (optional)
        line.erase(0, line.find_first_not_of(" \t\r\n")); // Trim leading whitespace
        line.erase(line.find_last_not_of(" \t\r\n") + 1); // Trim trailing whitespace
        if (line.empty())
        {
            continue;
        }

        size_t commentpos = line.find('#');
        if (commentpos > 0)
        {
            if (commentpos < std::string::npos)
            {
                line.erase(commentpos);
            }
        }
        else
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
    if (frame_start > frame_end-1)
    {
        std::cerr << "Error: Invalid frame range!" << std::endl;
        return;
    }
    
    // Load frame rate
    line_id ++;
    fps = std::stoi(lines[line_id]);

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
    int n_obj_class = 0;
    while (std::getline(parsed, line, ','))
    {
        line_id ++;
        if (line == "Tracer")
        {
            stb_list.push_back(STB<Tracer3D>(frame_start, frame_end, fps, vx_to_mm, n_thread, output_folder+"Tracer_"+std::to_string(n_obj_class)+'/', cam_list, axis_limit, lines[line_id]));

            // Calibrate OTF // 
            std::cout << "Start Calibrating OTF!" << std::endl;
            int n_otf_calib = 5;
            int n_obj2d_max = 1000;
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
                stb_list[n_obj_class]
            );

            std::cout << "Finish Calibrating OTF!\n" << std::endl;

            n_obj_class ++;
        }
        else
        {
            std::cerr << "Error: Unknown object type: " << line << std::endl;
            return;
        }
    }
    parsed.clear();

    // Load previous tracks
    line_id ++;
    int frame_id_prev = frame_start-1;

    if (line_id < lines.size())
    {
        std::cout << lines[line_id] << std::endl;
        bool is_load_tracks = false;

        parsed.str(lines[line_id]);
        std::getline(parsed, line, ',');
        is_load_tracks = std::stoi(line);
        std::getline(parsed, line, ',');
        frame_id_prev = std::stoi(line);
        parsed.clear();

        if (is_load_tracks)
        {
            if (frame_id_prev < frame_start-1 || frame_id_prev > frame_end-1)
            {
                std::cerr << "Error: Invalid previous frame ID!" << std::endl;
                return;
            }

            for (int i = 0; i < n_obj_class; i ++)
            {
                line_id ++;
                std::visit(
                    [&](auto& stb) 
                    { 
                        stb.loadTracks(lines[line_id], LONG_ACTIVE); 
                        stb.loadTracks(lines[line_id + n_obj_class], SHORT_ACTIVE);
                    }, 
                    stb_list[i]
                );
            }
        }
        else 
        {
            frame_id_prev = frame_start-1;
        }
    }

    // Start processing
    std::vector<Image> img_list(n_cam_all);
    for (int frame_id = frame_id_prev+1; frame_id < frame_end+1; frame_id ++)
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
    std::cout << "***************" << std::endl;
    std::cout << "OpenLPT finish!" << std::endl;
    std::cout << "***************" << std::endl;
}

int main (int argc, char *argv[])
{
    if (argc != 2)
    {
        std::cerr << "Main Error: Invalid number of arguments!" << std::endl;
        return 0;
    }

    std::string file = argv[1];
    run(file);

    return 0;
}