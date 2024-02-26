#include "test.h"
#include "ObjectFinder.h"
#include "ImageIO.h"
#include <string>
#include <iomanip>
#include <iostream>

bool test_function_1 ()
{
    ImageIO imgio;
    imgio.loadImgPath("../test/inputs/test_ObjectFinder/", "test_function_1.txt");
    Image img = imgio.loadImg(0);

    ObjectFinder2D objfinder;
    std::vector<Tracer2D> tr2d_list;
    std::vector<double> properties = {255, 30};
    objfinder.findObject2D(tr2d_list, img, properties);

    std::ofstream outfile("../test/results/test_ObjectFinder/test_function_1.csv");
    outfile.setf(std::ios_base::scientific);
    outfile.precision(8);
    for (auto& tr2d : tr2d_list)
    {
        outfile << tr2d._pt_center[0] << "," << tr2d._pt_center[1] << std::endl;
    }
    outfile.close();

    // compare with solution
    std::ifstream infile1("../test/solutions/test_ObjectFinder/test_function_1.csv");
    std::ifstream infile2("../test/results/test_ObjectFinder/test_function_1.csv");
    std::string line1, line2;
    while (std::getline(infile1, line1) && std::getline(infile2, line2))
    {
        if (line1 != line2)
        {
            std::cerr << "test_function_1() failed: tr2d_list is not correct!" << std::endl;
            std::cerr << "Expect:" << line1 << std::endl;
            std::cerr << "Get:" << line2 << std::endl;
            return false;
        }
    }
    infile1.close();
    infile2.close();

    return true;
}   

int main ()
{
    IS_TRUE(test_function_1());

    return 0;
}