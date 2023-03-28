#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <typeinfo>

#include <vector>
// #include <deque>
#include <typeinfo>
#include <omp.h>
#include <ctime>

#include "myMATH.hpp"
#include "Matrix.h"
#include "myIO.hpp"
// #include "Camera.h"
// #include "ImageIO.h"
// #include "ObjectFinder.h"
// #include "StereoMatch.h"
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

    std::string file_zsj = "Result/zsj.csv";
    std::string file_cyj = "Result/cyj.csv";

    // std::vector<std::vector<double>> zsj;
    // std::vector<std::vector<double>> cyj;
    // myIO::LoadMatrix<double> (file_zsj, zsj);
    // myIO::LoadMatrix<double> (file_cyj, cyj);
    Matrix<double> zsj(file_zsj);
    Matrix<double> cyj(file_cyj);

    Matrix<double> res = zsj * cyj;
    zsj *= cyj;

    zsj.WriteMatrix(std::string("Result/zsj_out.csv"));
    res.WriteMatrix(std::string("Result/cyj_out.csv"));
    // cyj.WriteMatrix(std::string("Result/cyj_out.csv"));

    std::cout << "Test ended!" << std::endl;
    return 0;
}
