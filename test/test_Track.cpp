#include "test.h"

#include <fstream>
#include <iostream>
#include <string>

#include "Track.h"
#include "Matrix.h"
#include "ObjectInfo.h"

bool test_function_1 ()
{
    Track<Tracer3D> track;

    Tracer3D obj3d;
    obj3d._pt_center[0] = 1;
    obj3d._pt_center[1] = 2;
    obj3d._pt_center[2] = 3;
    obj3d._error = 4;
    obj3d._n_2d = 1;
    obj3d._camid_list.push_back(6);
    obj3d._tr2d_list.push_back(Tracer2D(Pt2D(7, 8)));

    track.addNext(obj3d, 5);

    int n_cam_all = 10;
    int fps = 1;
    std::ofstream output("../test/results/test_Track/test_function_1.csv");
    output << "track_id,t,world_x,world_y,world_z,error";
    for (int i = 0; i < n_cam_all; i ++)
    {
        output << ",cam" << i << "_x(col),cam" << i << "_y(row)";
    }
    output << "\n";

    track.saveTrack(output, 0, fps, n_cam_all);
    track.saveTrack(output, 1, fps, n_cam_all);

    output.close();

    return true;
}


int main ()
{
    fs::create_directory("../test/results/test_Track/");

    IS_TRUE(test_function_1());

    return 0;
}