#include "test.h"
#include "OTF.h"


bool test_function_1()
{
    AxisLimit boundary;
    boundary.x_min = 0;
    boundary.x_max = 1;
    boundary.y_min = 0;
    boundary.y_max = 1;
    boundary.z_min = 0;
    boundary.z_max = 1;

    std::vector<double> otf_param_sol = {125,1.5,1.5,0};
    
    Pt3D pt3d(-10,0.5,10);
    OTF otf;
    otf.loadParam(1, 2, 3, 4, boundary);
    std::vector<double> otf_param = otf.getOTFParam(0, pt3d);

    if (otf_param != otf_param_sol)
    {
        std::cout << "test_function_1 failed" << std::endl;

        std::cout << "otf_param = ";
        for (int i = 0; i < otf_param.size(); i++)
        {
            std::cout << otf_param[i] << ',';
        }
        std::cout << std::endl;

        std::cout << "otf_param_sol = ";
        for (int i = 0; i < otf_param_sol.size(); i++)
        {
            std::cout << otf_param_sol[i] << ',';
        }
        std::cout << std::endl;

        return false;
    }

    return true;
}

bool test_function_2()
{
    std::vector<double> otf_param_sol = {125,1.5,1.5,0};
    
    Pt3D pt3d(-10,0.5,10);
    OTF otf;
    otf.loadParam("../test/inputs/test_OTF/otf.txt");
    std::vector<double> otf_param = otf.getOTFParam(0, pt3d);

    if (otf_param != otf_param_sol)
    {
        std::cout << "test_function_2 failed" << std::endl;

        std::cout << "otf_param = ";
        for (int i = 0; i < otf_param.size(); i++)
        {
            std::cout << otf_param[i] << ',';
        }
        std::cout << std::endl;

        std::cout << "otf_param_sol = ";
        for (int i = 0; i < otf_param_sol.size(); i++)
        {
            std::cout << otf_param_sol[i] << ',';
        }
        std::cout << std::endl;

        std::cout << "n_cam = " << otf._param.n_cam << std::endl;
        std::cout << "nx = " << otf._param.nx << std::endl;
        std::cout << "ny = " << otf._param.ny << std::endl;
        std::cout << "nz = " << otf._param.nz << std::endl;
        std::cout << "n_grid = " << otf._param.n_grid << std::endl;

        otf._param.a.print();
        otf._param.b.print();
        otf._param.c.print();
        otf._param.alpha.print();

        return false;
    }
    
    return true;
}

int main ()
{
    IS_TRUE(test_function_1());
    IS_TRUE(test_function_2());

    return 0;
}



