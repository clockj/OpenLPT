#include "test.h"

#include "Matrix.h"
#include "myMATH.h"


// test initialization
bool test_function_1 ()
{
    Matrix<double> ans("../test/solutions/test_Matrix/test_function_1.csv");
    std::vector<std::vector<double>> ans_vec = {{1,1,1},
                                                {1,1,1},
                                                {1,1,1},
                                                {1,2,3},
                                                {4,5,6},
                                                {7,8,9},
                                                {1,2,3},
                                                {4,5,6},
                                                {7,8,9}};
    for (int i = 0; i < ans.getDimX(); i ++)
    {
        for (int j = 0; j < ans.getDimY(); j ++)
        {
            if (ans(i,j) != ans_vec[i][j])
            {
                std::cout << "ans(" << i << "," << j << ") = " << ans(i,j) << ", ans = " << 0 << std::endl;
                return false;
            }
        }
    }
    int row_id = 0;

    Matrix<double> mtx_1(3, 3, 1);
    for (int i = 0; i < mtx_1.getDimX(); i ++)
    {
        for (int j = 0; j < mtx_1.getDimY(); j ++)
        {
            if (mtx_1(i,j) != ans(row_id,j))
            {
                std::cout << "mtx_1(" << i << "," << j << ") = " << mtx_1(i,j) << ", ans = " << ans(row_id,0) << std::endl;
                return false;
            }
        }
        row_id ++;
    }

    Matrix<double> mtx_2(3, 3, {{1,2,3},{4,5,6},{7,8,9}});
    for (int i = 0; i < mtx_2.getDimX(); i ++)
    {
        for (int j = 0; j < mtx_2.getDimY(); j ++)
        {
            if (mtx_2(i,j) != ans(row_id,j))
            {
                std::cout << "mtx_2(" << i << "," << j << ") = " << mtx_2(i,j) << ", ans = " << ans(row_id,0) << std::endl;
                return false;
            }
        }
        row_id ++;
    }

    Matrix<double> mtx_3(mtx_2);
    for (int i = 0; i < mtx_3.getDimX(); i ++)
    {
        for (int j = 0; j < mtx_3.getDimY(); j ++)
        {
            if (mtx_3(i,j) != ans(row_id,j))
            {
                std::cout << "mtx_3(" << i << "," << j << ") = " << mtx_3(i,j) << ", ans = " << ans(row_id,0) << std::endl;
                return false;
            }
        }
        row_id ++;
    }

    // test for large matrix
    Matrix<int> mtx_4(1024, 1024, 255);
    for (int i = 0; i < mtx_4.getDimX(); i ++)
    {
        for (int j = 0; j < mtx_4.getDimY(); j ++)
        {
            if (mtx_4(i,j) != 255)
            {
                std::cout << "mtx_4(" << i << "," << j << ") = " << mtx_4(i,j) << ", ans = " << 255 << std::endl;
                return false;
            }
        }
    }

    return true;
}

// test assign and get value 
bool test_function_2 ()
{
    Matrix<double> mtx_1(3,3, {{1,2,3},{4,5,6},{7,8,9}});
    if (mtx_1(1,1) != 5)
    {
        std::cout << "mtx_1(1,1) = " << mtx_1(1,1) << ", ans = " << 5 << std::endl;
        return false;
    }

    if (mtx_1[4] != 5)
    {
        std::cout << "mtx_1[4] = " << mtx_1[4] << ", ans = " << 5 << std::endl;
        return false;
    }

    mtx_1(1,1) = 100;
    if (mtx_1(1,1) != 100)
    {
        std::cout << "mtx_1(1,1) = " << mtx_1(1,1) << ", ans = " << 100 << std::endl;
        return false;
    }

    mtx_1[4] = 100;
    if (mtx_1[4] != 100)
    {
        std::cout << "mtx_1[4] = " << mtx_1[4] << ", ans = " << 100 << std::endl;
        return false;
    }

    return true;
}

// test get row and col
bool test_function_3 ()
{
    Matrix<double> mtx_1(3,3, {{1,2,3},{4,5,6},{7,8,9}});

    std::vector<double> row = mtx_1.getRow(1);
    std::vector<double> ans = {4,5,6};
    if (row.size() != ans.size())
    {
        std::cout << "row.size = " << row.size() << ", ans = " << ans.size() << std::endl;
        return false;
    }
    for (int i = 0; i < row.size(); i ++)
    {
        if (row[i] != ans[i])
        {
            std::cout << "row[" << i << "] = " << row[i] << ", ans = " << ans[i] << std::endl;
            return false;
        }
    }

    std::vector<double> col = mtx_1.getCol(1);
    ans = {2,5,8};
    if (col.size() != ans.size())
    {
        std::cout << "col.size = " << col.size() << ", ans = " << ans.size() << std::endl;
        return false;
    }
    for (int i = 0; i < col.size(); i ++)
    {
        if (col[i] != ans[i])
        {
            std::cout << "col[" << i << "] = " << col[i] << ", ans = " << ans[i] << std::endl;
            return false;
        }
    }

    return true;
}

// test get matrix info
bool test_function_4 ()
{
    Matrix<double> mtx_1(100,500, 10);
    if (mtx_1.getDimX() != 100)
    {
        std::cout << "mtx_1.getDimX() = " << mtx_1.getDimX() << ", ans = " << 100 << std::endl;
        return false;
    }
    if (mtx_1.getDimY() != 500)
    {
        std::cout << "mtx_1.getDimY() = " << mtx_1.getDimY() << ", ans = " << 500 << std::endl;
        return false;
    }

    return true;
}

// test magnitude operations
bool test_function_5 ()
{
    Matrix<double> mtx_1(3,3, {{1,2,3},{4,5,6},{7,8,9}});
    double norm = mtx_1.norm();
    std::stringstream ss;
    ss << std::fixed << std::setprecision(4) << norm;
    ss >> norm;
    if (norm != 16.8819)
    {
        std::cout << "mtx_1.norm() = " << norm << ", ans = " << 16.8819 << std::endl;
        return false;
    }
    
    return true;
}

// test matrix operations 
bool test_function_6 ()
{
    Matrix<double> mtx_1(3,3, {{1,2,3},{4,5,6},{7,8,9}});

    Matrix<double> mtx_2 = mtx_1;
    for (int i = 0; i < mtx_2.getDimX(); i ++)
    {
        for (int j = 0; j < mtx_2.getDimY(); j ++)
        {
            if (mtx_2(i,j) != mtx_1(i,j))
            {
                std::cout << "mtx_2(" << i << "," << j << ") = " << mtx_2(i,j) << ", ans = " << mtx_1(i,j) << std::endl;
                return false;
            }
        }
    }

    Matrix<double> mtx_3 = mtx_1 + mtx_2;
    for (int i = 0; i < mtx_3.getDimX(); i ++)
    {
        for (int j = 0; j < mtx_3.getDimY(); j ++)
        {
            if (mtx_3(i,j) != 2*mtx_1(i,j))
            {
                std::cout << "mtx_3(" << i << "," << j << ") = " << mtx_3(i,j) << ", ans = " << 2*mtx_1(i,j) << std::endl;
                return false;
            }
        }
    }

    Matrix<double> mtx_4 = mtx_1 - mtx_2;
    for (int i = 0; i < mtx_4.getDimX(); i ++)
    {
        for (int j = 0; j < mtx_4.getDimY(); j ++)
        {
            if (mtx_4(i,j) != 0)
            {
                std::cout << "mtx_4(" << i << "," << j << ") = " << mtx_4(i,j) << ", ans = " << 0 << std::endl;
                return false;
            }
        }
    }

    Matrix<double> mtx_5 = mtx_1 * mtx_2;
    Matrix<double> ans_5(3, 3, {{30,36,42},{66,81,96},{102,126,150}});
    for (int i = 0; i < mtx_5.getDimX(); i ++)
    {
        for (int j = 0; j < mtx_5.getDimY(); j ++)
        {
            if (mtx_5(i,j) != ans_5(i,j))
            {
                std::cout << "mtx_5(" << i << "," << j << ") = " << mtx_5(i,j) << ", ans = " << 30 << std::endl;
                return false;
            }
        }
    }

    Matrix<double> mtx_6 = mtx_1 * 2;
    for (int i = 0; i < mtx_6.getDimX(); i ++)
    {
        for (int j = 0; j < mtx_6.getDimY(); j ++)
        {
            if (mtx_6(i,j) != 2*mtx_1(i,j))
            {
                std::cout << "mtx_6(" << i << "," << j << ") = " << mtx_6(i,j) << ", ans = " << 2*mtx_1(i,j) << std::endl;
                return false;
            }
        }
    }

    Matrix<double> mtx_7 = mtx_1 / 2;
    for (int i = 0; i < mtx_7.getDimX(); i ++)
    {
        for (int j = 0; j < mtx_7.getDimY(); j ++)
        {
            if (mtx_7(i,j) != mtx_1(i,j)/2)
            {
                std::cout << "mtx_7(" << i << "," << j << ") = " << mtx_7(i,j) << ", ans = " << mtx_1(i,j)/2 << std::endl;
                return false;
            }
        }
    }

    Matrix<double> mtx_8 = mtx_1;
    mtx_8 += mtx_2;
    for (int i = 0; i < mtx_8.getDimX(); i ++)
    {
        for (int j = 0; j < mtx_8.getDimY(); j ++)
        {
            if (mtx_8(i,j) != mtx_1(i,j) + mtx_2(i,j))
            {
                std::cout << "mtx_8(" << i << "," << j << ") = " << mtx_8(i,j) << ", ans = " << 2*mtx_1(i,j) << std::endl;
                return false;
            }
        }
    }

    Matrix<double> mtx_9 = mtx_1;
    mtx_9 -= mtx_2;
    for (int i = 0; i < mtx_9.getDimX(); i ++)
    {
        for (int j = 0; j < mtx_9.getDimY(); j ++)
        {
            if (mtx_9(i,j) != 0)
            {
                std::cout << "mtx_9(" << i << "," << j << ") = " << mtx_9(i,j) << ", ans = " << 0 << std::endl;
                return false;
            }
        }
    }

    Matrix<double> mtx_10 = mtx_1;
    mtx_10 *= mtx_2;
    for (int i = 0; i < mtx_10.getDimX(); i ++)
    {
        for (int j = 0; j < mtx_10.getDimY(); j ++)
        {
            if (mtx_10(i,j) != ans_5(i,j))
            {
                std::cout << "mtx_10(" << i << "," << j << ") = " << mtx_10(i,j) << ", ans = " << 30 << std::endl;
                return false;
            }
        }
    }

    Matrix<double> mtx_11 = mtx_1;
    mtx_11 *= 2;
    for (int i = 0; i < mtx_11.getDimX(); i ++)
    {
        for (int j = 0; j < mtx_11.getDimY(); j ++)
        {
            if (mtx_11(i,j) != 2*mtx_1(i,j))
            {
                std::cout << "mtx_11(" << i << "," << j << ") = " << mtx_11(i,j) << ", ans = " << 2*mtx_1(i,j) << std::endl;
                return false;
            }
        }
    }

    Matrix<double> mtx_12 = mtx_1;
    mtx_12 /= 2;
    for (int i = 0; i < mtx_12.getDimX(); i ++)
    {
        for (int j = 0; j < mtx_12.getDimY(); j ++)
        {
            if (mtx_12(i,j) != mtx_1(i,j)/2)
            {
                std::cout << "mtx_12(" << i << "," << j << ") = " << mtx_12(i,j) << ", ans = " << mtx_1(i,j)/2 << std::endl;
                return false;
            }
        }
    }

    Matrix<double> mtx_13(3,4, {{1,2,3,4},{5,6,7,8},{9,10,11,12}});
    Matrix<double> mtx_14(4,3, {{1,2,3},{4,5,6},{7,8,9},{10,11,12}});
    Matrix<double> mtx_15 = mtx_13 * mtx_14;
    Matrix<double> ans_15(3,3, {{70,80,90},{158,184,210},{246,288,330}});
    for (int i = 0; i < mtx_15.getDimX(); i ++)
    {
        for (int j = 0; j < mtx_15.getDimY(); j ++)
        {
            if (mtx_15(i,j) != ans_15(i,j))
            {
                std::cout << "mtx_15(" << i << "," << j << ") = " << mtx_15(i,j) << ", ans = " << ans_15(i,j) << std::endl;
                return false;
            }
        }
    }

    Matrix<double> mtx_16(3,3, {{1,2,3},{4,5,6},{7,8,9}});
    Matrix<double> mtx_17 = mtx_16.transpose();
    Matrix<double> ans_17(3,3, {{1,4,7},{2,5,8},{3,6,9}});
    for (int i = 0; i < mtx_17.getDimX(); i ++)
    {
        for (int j = 0; j < mtx_17.getDimY(); j ++)
        {
            if (mtx_17(i,j) != ans_17(i,j))
            {
                std::cout << "mtx_17(" << i << "," << j << ") = " << mtx_17(i,j) << ", ans = " << ans_17(i,j) << std::endl;
                return false;
            }
        }
    }

    Matrix<double> mtx_18(3,3, {{1,2,3},{4,5,6},{7,8,9}});
    Matrix<double> mtx_19 = mtx_18;
    if (mtx_18 != mtx_19)
    {
        std::cout << "result: mtx_18 != mtx_19, answer: mtx_18 == mtx_19" << std::endl;
        return false;
    }
    if (!(mtx_18 == mtx_1))
    {
        std::cout << "result: mtx_18 != mtx_1, answer: mtx_18 == mtx_1" << std::endl;
        return false;
    }
    if (mtx_18 != mtx_18)
    {
        std::cout << "result: mtx_18 != mtx_18, answer: mtx_18 == mtx_18" << std::endl;
        return false;
    }

    return true;
}

// test matrix hierarchy
bool test_function_7 ()
{
    Pt3D pt_1(1,2,3);
    for (int i = 0; i < 3; i ++)
    {
        if (pt_1[i] != i+1)
        {
            std::cout << "pt_1[" << i << "] = " << pt_1[i] << ", ans = " << i+1 << std::endl;
            return false;
        }
    }

    pt_1 *= 3;
    for (int i = 0; i < 3; i ++)
    {
        if (pt_1[i] != 3*(i+1))
        {
            std::cout << "pt_1[" << i << "] = " << pt_1[i] << ", ans = " << 3*(i+1) << std::endl;
            return false;
        }
    }

    Matrix<double> mtx_1(3,3, {{1,2,3},{4,5,6},{7,8,9}});
    Pt3D pt_2 = (mtx_1 * pt_1);
    std::vector<double> ans = {42,96,150};
    for (int i = 0; i < 3; i ++)
    {
        if (pt_2[i] != ans[i])
        {
            std::cout << "pt_2[" << i << "] = " << pt_2[i] << ", ans = " << ans[i] << std::endl;
            return false;
        }
    }

    Pt2D pt_3(1,2);
    for (int i = 0; i < 2; i ++)
    {
        if (pt_3[i] != i+1)
        {
            std::cout << "pt_3[" << i << "] = " << pt_3[i] << ", ans = " << i+1 << std::endl;
            return false;
        }
    }

    pt_3 *= 3;
    for (int i = 0; i < 2; i ++)
    {
        if (pt_3[i] != 3*(i+1))
        {
            std::cout << "pt_3[" << i << "] = " << pt_3[i] << ", ans = " << 3*(i+1) << std::endl;
            return false;
        }
    }

    Matrix<double> mtx_2(2,2, {{1,2},{3,4}});
    Pt2D pt_4 = (mtx_2 * pt_3);
    ans = {15,33};
    for (int i = 0; i < 2; i ++)
    {
        if (pt_4[i] != ans[i])
        {
            std::cout << "pt_4[" << i << "] = " << pt_4[i] << ", ans = " << ans[i] << std::endl;
            return false;
        }
    }

    Image img_1(1024,1024, 255);
    for (int i = 0; i < img_1.getDimX(); i ++)
    {
        for (int j = 0; j < img_1.getDimY(); j ++)
        {
            if (img_1(i,j) != 255)
            {
                std::cout << "img_1(" << i << "," << j << ") = " << img_1(i,j) << ", ans = " << 255 << std::endl;
                return false;
            }
        }
    }

    Matrix<double> mtx_3(1024,1024, 1);
    img_1 = mtx_3;
    for (int i = 0; i < img_1.getDimX(); i ++)
    {
        for (int j = 0; j < img_1.getDimY(); j ++)
        {
            if (img_1(i,j) != 1)
            {
                std::cout << "img_1(" << i << "," << j << ") = " << img_1(i,j) << ", ans = " << 1 << std::endl;
                return false;
            }
        }
    }

    return true;
}

int main()
{
    IS_TRUE(test_function_1());
    IS_TRUE(test_function_2());
    IS_TRUE(test_function_3());
    IS_TRUE(test_function_4());
    IS_TRUE(test_function_5());
    IS_TRUE(test_function_6());
    IS_TRUE(test_function_7());
    
    return 0;
}