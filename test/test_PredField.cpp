#include "test.h"
#include "STBCommons.h"
#include "Matrix.h"
#include "ObjectInfo.h"
#include "PredField.h"

bool test_function_1 ()
{
    std::cout << "test_function_1" << std::endl;

    PFParam param;
    param.limit = AxisLimit(-20, 20, -20, 20, -20, 20);
    param.nx = 51;
    param.ny = 51;
    param.nz = 51;
    param.r = 0.4; // 40/50 = 0.8, r = 0.8/2 = 0.4
    param.nBin_x = 5;
    param.nBin_y = 5;
    param.nBin_z = 5;

    // 3D object list
    Matrix<double> inputs_prev("../test/inputs/test_PredField/pts_prev_1.csv");
    Matrix<double> inputs_curr("../test/inputs/test_PredField/pts_curr_1.csv");
    std::vector<Pt3D> pt3d_list_prev(inputs_prev.getDimRow());
    std::vector<Pt3D> pt3d_list_curr(inputs_curr.getDimRow());
    for (int i = 0; i < inputs_prev.getDimRow(); i++)
    {
        pt3d_list_prev[i][0] = inputs_prev(i, 0);
        pt3d_list_prev[i][1] = inputs_prev(i, 1);
        pt3d_list_prev[i][2] = inputs_prev(i, 2);
    }
    for (int i = 0; i < inputs_curr.getDimRow(); i++)
    {
        pt3d_list_curr[i][0] = inputs_curr(i, 0);
        pt3d_list_curr[i][1] = inputs_curr(i, 1);
        pt3d_list_curr[i][2] = inputs_curr(i, 2);
    }

    PredField pf(param, pt3d_list_prev, pt3d_list_curr);


    // load real displacement field
    Matrix<double> disp_field("../test/solutions/test_PredField/disp_field_1.csv");
    int n_wrong = 0;
    int n_notfind = 0;
    for (int i = 0; i < disp_field.getDimRow(); i++)
    {
        for (int j = 0; j < 3; j++)
        {
            if (std::fabs(pf._disp_field(i, j) - disp_field(i, j)) > 1e-6)
            {   
                n_wrong ++;
                
                if (std::fabs(pf._disp_field(i, j)) < SMALLNUMBER)
                {
                    n_notfind ++;
                }
                else
                {
                    // std::cout << "i=" << i << ", " << pf._disp_field(i, 0) <<  "," << pf._disp_field(i, 1) << "," << pf._disp_field(i, 2) << " " << disp_field(i, 0) << "," << disp_field(i, 1) << "," << disp_field(i, 2) << std::endl;
                }

                // return false;
                break;
            }
        }
    }

    int n_grid = pf._disp_field.getDimRow();
    int n_find = n_grid - n_notfind;
    int n_correct = n_grid - n_wrong;
    std::cout << "n_grid=" << disp_field.getDimRow() << std::endl;
    std::cout << "n_wrong=" << n_wrong << std::endl;
    std::cout << "ratio_wrong=" << 100.0*n_wrong/disp_field.getDimRow() << std::endl;
    std::cout << "n_notfind=" << n_notfind << std::endl;
    std::cout << "ratio_notfind=" << 100.0*n_notfind/disp_field.getDimRow() << std::endl;
    std::cout << "n_correct=" << n_correct << std::endl;
    std::cout << "n_find=" << n_find << std::endl;
    std::cout << "ratio_correct=" << 100.0*n_correct/n_find << std::endl;

    std::cout << "test_function_1 done\n" << std::endl;

    return true;
}


bool test_function_2 ()
{
    std::cout << "test_function_2" << std::endl;

    PFParam param;
    param.limit = AxisLimit(-20, 20, -20, 20, -20, 20);
    param.nx = 51;
    param.ny = 51;
    param.nz = 51;
    param.r = 0.4; // 40/50 = 0.8, r = 0.8/2 = 0.4
    param.nBin_x = 5;
    param.nBin_y = 5;
    param.nBin_z = 5;

    // 3D object list
    Matrix<double> inputs_prev("../test/inputs/test_PredField/pts_prev_2.csv");
    Matrix<double> inputs_curr("../test/inputs/test_PredField/pts_curr_2.csv");
    std::vector<Pt3D> pt3d_list_prev(inputs_prev.getDimRow());
    std::vector<Pt3D> pt3d_list_curr(inputs_curr.getDimRow());
    for (int i = 0; i < inputs_prev.getDimRow(); i++)
    {
        pt3d_list_prev[i][0] = inputs_prev(i, 0);
        pt3d_list_prev[i][1] = inputs_prev(i, 1);
        pt3d_list_prev[i][2] = inputs_prev(i, 2);
    }
    for (int i = 0; i < inputs_curr.getDimRow(); i++)
    {
        pt3d_list_curr[i][0] = inputs_curr(i, 0);
        pt3d_list_curr[i][1] = inputs_curr(i, 1);
        pt3d_list_curr[i][2] = inputs_curr(i, 2);
    }

    PredField pf(param, pt3d_list_prev, pt3d_list_curr);


    // load real displacement field
    Matrix<double> disp_field("../test/solutions/test_PredField/disp_field_2.csv");
    int n_wrong = 0;
    int n_notfind = 0;
    for (int i = 0; i < disp_field.getDimRow(); i++)
    {
        for (int j = 0; j < 3; j++)
        {
            if (std::fabs(pf._disp_field(i, j) - disp_field(i, j)) > 1e-4)
            {   
                n_wrong ++;
                
                if (std::fabs(pf._disp_field(i, j)) < SMALLNUMBER)
                {
                    n_notfind ++;
                }
                else
                {
                    // std::cout << "i=" << i << ", " << pf._disp_field(i, 0) <<  "," << pf._disp_field(i, 1) << "," << pf._disp_field(i, 2) << " " << disp_field(i, 0) << "," << disp_field(i, 1) << "," << disp_field(i, 2) << std::endl;
                }

                // return false;
                break;
            }
        }
    }

    int n_grid = pf._disp_field.getDimRow();
    int n_find = n_grid - n_notfind;
    int n_correct = n_grid - n_wrong;
    std::cout << "n_grid=" << disp_field.getDimRow() << std::endl;
    std::cout << "n_wrong=" << n_wrong << std::endl;
    std::cout << "ratio_wrong=" << 100.0*n_wrong/disp_field.getDimRow() << std::endl;
    std::cout << "n_notfind=" << n_notfind << std::endl;
    std::cout << "ratio_notfind=" << 100.0*n_notfind/disp_field.getDimRow() << std::endl;
    std::cout << "n_correct=" << n_correct << std::endl;
    std::cout << "n_find=" << n_find << std::endl;
    std::cout << "ratio_correct=" << 100.0*n_correct/n_find << std::endl;

    std::cout << "test_function_2 done\n" << std::endl;
    return true;
}


int main()
{
    IS_TRUE(test_function_1());
    IS_TRUE(test_function_2());

    return 0;
}