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
// #include "Matrix.h"
// #include "Camera.h"
// #include "ImageIO.h"
// #include "myMATH.h"
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

void WriteMatchList (std::string file_name, std::vector<std::vector<int>>& tracer_match_list)
{
    std::cout << "Start writing match list!" << std::endl;

    std::ofstream outfile(file_name, std::ios::out);

    outfile << "MatchList" << "\n";
    outfile << "Number of matches is: " << tracer_match_list.size() << "\n";

    for (int i = 0; i < tracer_match_list.size(); i ++)
    {
        for (int j = 0; j < tracer_match_list[i].size(); j ++)
        {
            outfile << tracer_match_list[i][j] << "\t";
        }
        outfile << "\n";
    }

    std::cout << "Finish writing match list!" << std::endl;
}

void WriteTriError (std::string file_name, std::vector<double>& error_list)
{
    std::cout << "Start writing tri-error list!" << std::endl;

    std::ofstream outfile(file_name, std::ios::out);

    outfile << "MatchList" << "\n";
    outfile << "Number of matches is: " << error_list.size() << "\n";

    for (int i = 0; i < error_list.size(); i ++)
    {
        outfile << error_list[i] << "\n";
    }

    std::cout << "Finish writing tri-error list!" << std::endl;
}


int main()
{
    std::cout << "Simple test!" << std::endl;

    std::vector<double> nums = {-0.9,0,0.5,0.5,0.1,800,0.32,900,900,9990,1231};
    std::vector<int> sorted_index(nums.size());

    myMATH::MergeSort<double>(sorted_index, nums);
    // myMATH::BubbleSort<double>(sorted_index, nums);
    for (int i = 0; i < nums.size(); i ++)
    {
        std::cout << nums[sorted_index[i]] << " ";
    }
    std::cout << std::endl;

    std::cout << "Max=" << myMATH::Max<double>(nums) << std::endl;
    std::cout << "Min=" << myMATH::Min<double>(nums) << std::endl;
    std::cout << "Medium=" << myMATH::Median<double>(nums) << std::endl;

    std::vector<bool> judge(nums.size());
    myMATH::IsOutlier(judge, nums);
    for (int i = 0; i < nums.size(); i ++)
    {
        std::cout << judge[i] << " ";
    }
    std::cout << std::endl;

    std::cout << error_type << std::endl;

    std::cout << "Test ended!" << std::endl;
    return 0;
}
