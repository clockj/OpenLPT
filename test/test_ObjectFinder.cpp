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
    std::vector<double> properties = {255, 30, 2};
    objfinder.findObject2D(tr2d_list, img, properties);

    std::ofstream outfile("../test/results/test_ObjectFinder/test_function_1.csv");
    outfile.setf(std::ios_base::scientific);
    outfile.precision(SAVEPRECISION);
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

// compare with real 2d location
bool test_function_2 ()
{
    // load image path
    std::vector<ImageIO> imgio_list;
    for (int i = 0; i < 4; i ++)
    {
        ImageIO imgio;
        imgio.loadImgPath("../test/inputs/test_ObjectFinder/", "cam" + std::to_string(i+1) + "ImageNames" + ".txt");
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
    std::vector<double> properties = {255, 30};
    ObjectFinder2D objfinder;
    for (int i = 0; i < 4; i ++)
    {
        std::vector<Tracer2D> tr2d_list;
        objfinder.findObject2D(tr2d_list, img_list[i], properties);
        tr2d_list_all.push_back(tr2d_list);
        std::cout << "tr2d_list_all[" << i << "].size() = " << tr2d_list_all[i].size() << std::endl;

        // save result
        std::ofstream outfile("../test/results/test_ObjectFinder/pt2d_list_cam" + std::to_string(i+1) + ".csv");
        outfile.precision(SAVEPRECISION);
        for (auto& tr2d : tr2d_list)
        {
            outfile << tr2d._pt_center[0] << "," << tr2d._pt_center[1] << std::endl;
        }
    }

    return true;
}

int main ()
{
    IS_TRUE(test_function_1());
    IS_TRUE(test_function_2());

    return 0;
}