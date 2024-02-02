#include <vector>
// #include <deque>
#include <typeinfo>
#include <omp.h>
#include <ctime>
#include <cstring>

#include "myMATH.h"
#include "Matrix.h"
// #include "Camera.h"
// #include "ObjectInfo.h"
#include "ImageIO.h"
// #include "ObjectFinder.h"
// #include "StereoMatch.h"
// #include "OTF.h"
// #include "Shake.h"
// #include "IPR.h"
// #include "PredField.h"
// #include "Track.h"
// #include "STB.h"

void LoadPtList(std::vector<Matrix<double>>& pt_list, std::string file_name)
{
    std::string line;
    int n_pt;
    double value;

    std::ifstream infile;
    std::istringstream istream;
    infile.open(file_name);

    std::getline(infile, line);
    istream.str(line);
    istream >> n_pt;

    pt_list = std::vector<Matrix<double>>(n_pt, Matrix<double>(3,1,0));

    for (int i = 0; i < n_pt; i ++)
    {
        std::getline(infile, line);
        std::istringstream istream_data;
        istream_data.str(line);
        for (int j = 0; j < 3; j ++)
        {   
            istream_data >> value;
            pt_list[i](j,0) = value;
            if (istream_data.peek() == ',')
            {
                istream_data.ignore();
            }
        }

    }
    infile.close();
};


int main()
{
    std::cout << "Simple test!" << std::endl;

    // // for debug
    // // OTF 
    // AxisLimit limit;
    // limit.x_min = -20; limit.x_max = 20;
    // limit.y_min = -20; limit.y_max = 20;
    // limit.z_min = -20; limit.z_max = 20;
    // OTF otf(4, 3, limit);

    // // Predicted field
    // std::cout << "PredField start!" << std::endl;
    // std::vector<std::vector<Matrix<double>>> pt_list_all;
    // // Matrix<double> mtx_pre("3Dpoints_pf_1.csv");
    // Matrix<double> mtx_pre("3Dpoints_ipr_1.csv");
    // std::vector<Matrix<double>> pt_list_pre(mtx_pre.GetDimX(), Matrix<double>(3,1,0));
    // for (int i = 0; i < mtx_pre.GetDimX(); i ++)
    // {
    //     pt_list_pre[i](0,0) = mtx_pre(i,0);
    //     pt_list_pre[i](1,0) = mtx_pre(i,1);
    //     pt_list_pre[i](2,0) = mtx_pre(i,2);
    // }
    // std::cout << pt_list_pre.size() << std::endl;
    // // Matrix<double> mtx_cur("3Dpoints_pf_2.csv");
    // Matrix<double> mtx_cur("3Dpoints_ipr_2.csv");
    // std::vector<Matrix<double>> pt_list_cur(mtx_cur.GetDimX(), Matrix<double>(3,1,0));
    // for (int i = 0; i < mtx_cur.GetDimX(); i ++)
    // {
    //     pt_list_cur[i](0,0) = mtx_cur(i,0);
    //     pt_list_cur[i](1,0) = mtx_cur(i,1);
    //     pt_list_cur[i](2,0) = mtx_cur(i,2);
    // }
    // std::cout << pt_list_cur.size() << std::endl;
    // pt_list_all.push_back(pt_list_pre);
    // pt_list_all.push_back(pt_list_cur);

    // std::vector<int> n_xyz = {21,21,21};
    // double r = 25 * 40/1000; // 25 * (max-min)/1000
    // std::vector<TracerInfo> useless;
    // PredField<TracerInfo> pf (limit, n_xyz, pt_list_all[0], pt_list_all[1], r, useless);
    // pf.SaveField("pf_disp_field.csv");

    // Matrix<double> disp_pred(3,1);
    // disp_pred = pf.PtInterp(pt_list_all[1][0]);
    // disp_pred.Print();

    // STB<TracerInfo> stb("./Data/stbConfig.txt");
    // // STB<TracerInfo> stb("D:/SD00125_New/stbConfig.txt");
    // // STB<TracerInfo> stb("D:/SD00125_New_updatesearch/stbConfig.txt");
    // stb.Run();

    // //Holo data analysis
    // int n_frame = 100;
    // std::vector<std::vector<Matrix<double>>> pt_list_all(n_frame);
    // for (int i = 0; i < n_frame; i ++)
    // {
    //     // LoadPtList(pt_list_all[i], "D:/My Code/Tracking Code/Holography track/Case/Data/"+std::to_string(i+1)+".csv");
    //     LoadPtList(pt_list_all[i], "D:/My Code/Tracking Code/Holography track/Case/Data_ml/"+std::to_string(i+1)+".csv");
    // }
    // std::cout << "Finish loading!" << std::endl;
    // std::vector<ObjectInfo> useless;
    // AxisLimit limit;
    // limit._x_min = 1;
    // limit._x_max = 306;
    // limit._y_min = 1;
    // limit._y_max = 882;
    // limit._z_min = 1;
    // limit._z_max = 1001;
    // std::vector<int> n_xyz = {10,28,32};
    // for (int i = 0; i < n_frame-1; i ++)
    // {
    //     PredField<ObjectInfo> pf(limit, n_xyz, pt_list_all[i], pt_list_all[i+1], 10, useless);
    //     // pf.SaveField("D:/My Code/Tracking Code/Holography track/Case/DispField/"+std::to_string(i+1)+".csv");
    //     pf.SaveField("D:/My Code/Tracking Code/Holography track/Case/DispField_ml/"+std::to_string(i+1)+".csv");
    // }

    std::cout << "Finish OpenLPT!" << std::endl;
    return 0;
}
