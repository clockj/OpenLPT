#include "test.h"
#include "myMATH.h"

// test sortID
bool test_function_1 ()
{   
    std::vector<double> nums = {5,4,1,2,3};
    std::vector<int> sort_index(nums.size());
    myMATH::sortID<double> (sort_index, nums);

    std::vector<int> sort_index_ans = {2,3,4,1,0};
    if (sort_index != sort_index_ans)
    {
        std::cout << "test_function_1: sortID failed" << std::endl;
        std::cout << "sort_index: " << std::endl;
        for (int i = 0; i < sort_index.size(); i ++)
        {
            std::cout << sort_index[i] << " ";
        }
        std::cout << std::endl;
        std::cout << "sort_index_ans: " << std::endl;
        for (int i = 0; i < sort_index_ans.size(); i ++)
        {
            std::cout << sort_index_ans[i] << " ";
        }
        std::cout << std::endl;
        return false;
    }

    return true;
}

// test getMedian
bool test_function_2 ()
{
    std::vector<double> nums = {5,4,1,2,3};
    double median = myMATH::getMedian<double> (nums);

    if (median != 3)
    {
        std::cout << "test_function_2: getMedian failed" << std::endl;
        std::cout << "median: " << median << std::endl;
        return false;
    }

    return true;
}

// test isOutlier
bool test_function_3 ()
{
    std::vector<double> nums = {10,4,1,2,3};
    std::vector<int> judge(nums.size());
    myMATH::isOutlier<double> (judge, nums);

    std::vector<int> judge_ans = {1,0,0,0,0};
    if (judge != judge_ans)
    {
        std::cout << "test_function_3: isOutlier failed" << std::endl;
        std::cout << "judge: " << std::endl;
        for (int i = 0; i < judge.size(); i ++)
        {
            std::cout << judge[i] << " ";
        }
        std::cout << std::endl;
        std::cout << "judge_ans: " << std::endl;
        for (int i = 0; i < judge_ans.size(); i ++)
        {
            std::cout << judge_ans[i] << " ";
        }
        std::cout << std::endl;
        return false;
    }

    return true;
}

// test linspace
bool test_function_4 ()
{
    std::vector<double> nums = myMATH::linspace (0, 1, 11);

    std::vector<double> nums_ans = {0,0.1,0.2,0.3,0.4,0.5,
                                    0.6,0.7,0.8,0.9,1};

    for (int i = 0; i < nums.size(); i ++)
    {
        if (std::fabs(nums[i] - nums_ans[i]) > SMALLNUMBER)
        {
            std::cout << "test_function_4: linspace failed" << std::endl;
            std::cout << "nums: " << std::endl;
            for (int i = 0; i < nums.size(); i ++)
            {
                std::cout << nums[i] << " ";
            }
            std::cout << std::endl;
            std::cout << "nums_ans: " << std::endl;
            for (int i = 0; i < nums_ans.size(); i ++)
            {
                std::cout << nums_ans[i] << " ";
            }
            std::cout << std::endl;
            return false;
        }
    }

    return true;
}

// test triLinearInterp
bool test_function_5 ()
{
    AxisLimit grid_limit;

    grid_limit.x_min = 0;
    grid_limit.x_max = 1;
    grid_limit.y_min = 0;
    grid_limit.y_max = 1;
    grid_limit.z_min = 0;
    grid_limit.z_max = 1;

    std::vector<double> value = {0,1,2,3,4,5,6,7};
    
    Matrix<double> inputs ("../test/inputs/test_myMATH/test_function_5.csv");
    Matrix<double> results (inputs.getDimRow(), 1, 0);
    for (int i = 0; i < inputs.getDimRow(); i ++)
    {
        std::vector<double> pt_vec = inputs.getRow(i);
        results[i] = myMATH::triLinearInterp(grid_limit, value, pt_vec);
    }
    results.write("../test/results/test_myMATH/test_function_5.csv");

    Matrix<double> answers ("../test/solutions/test_myMATH/test_function_5.csv");

    if (results != answers)
    {
        std::cout << "test_function_5: triLinearInterp failed" << std::endl;
        std::cout << "results: " << std::endl;
        results.print();
        std::cout << "answers: " << std::endl;
        answers.print();
        return false;
    }

    return true;
}

// test createUnitVector
bool test_function_6 ()
{
    Pt3D pt1(1,0,0);
    Pt3D pt2(0,1,0);
    Pt3D unit_vector = myMATH::createUnitVector (pt1, pt2);
    Pt3D unit_vector_ans(-1/std::sqrt(2),1/std::sqrt(2),0);

    if (unit_vector != unit_vector_ans)
    {
        std::cout << "test_function_6: createUnitVector failed" << std::endl;
        std::cout << "unit_vector 3D: " << std::endl;
        unit_vector.print();
        std::cout << "unit_vector_ans 3D: " << std::endl;
        unit_vector_ans.print();
        return false;
    }

    Pt2D pt3(1,0);
    Pt2D pt4(0,1);
    Pt2D unit_vector_2D = myMATH::createUnitVector (pt3, pt4);
    Pt2D unit_vector_ans_2D(-1/std::sqrt(2),1/std::sqrt(2));

    if (unit_vector_2D != unit_vector_ans_2D)
    {
        std::cout << "test_function_6: createUnitVector failed" << std::endl;
        std::cout << "unit_vector 2D: " << std::endl;
        unit_vector_2D.print();
        std::cout << "unit_vector_ans 2D: " << std::endl;
        unit_vector_ans_2D.print();
        return false;
    }

    return true;
}

// test dot product
bool test_function_7 ()
{
    Pt3D pt1(1,2,3);
    Pt3D pt2(4,5,6);
    double dot_product = myMATH::dot (pt1, pt2);
    double dot_product_ans = 32;

    if (std::fabs(dot_product - dot_product_ans) > SMALLNUMBER)
    {
        std::cout << "test_function_7: dot product 3D failed" << std::endl;
        std::cout << "dot_product 3D: " << dot_product << std::endl;
        std::cout << "dot_product_ans 3D: " << dot_product_ans << std::endl;
        return false;
    }

    Pt2D pt3(1,2);
    Pt2D pt4(4,5);
    double dot_product_2D = myMATH::dot (pt3, pt4);
    double dot_product_ans_2D = 14;

    if (std::fabs(dot_product_2D - dot_product_ans_2D) > SMALLNUMBER)
    {
        std::cout << "test_function_7: dot product 2D failed" << std::endl;
        std::cout << "dot_product_2D: " << dot_product_2D << std::endl;
        std::cout << "dot_product_ans_2D: " << dot_product_ans_2D << std::endl;
        return false;
    }

    return true;
}

// test calculate the distance between two points 
bool test_function_8 ()
{
    Pt3D pt1(1,2,3);
    Pt3D pt2(4,5,6);
    double distance_3D = myMATH::dist (pt1, pt2);
    double distance_3D_ans = std::sqrt(27);

    if (std::fabs(distance_3D - distance_3D_ans) > SMALLNUMBER)
    {
        std::cout << "test_function_8: distance 3D failed" << std::endl;
        std::cout << "distance_3D: " << distance_3D << std::endl;
        std::cout << "distance_3D_ans: " << distance_3D_ans << std::endl;
        return false;
    }

    Pt2D pt3(1,2);
    Pt2D pt4(4,5);
    double distance_2D = myMATH::dist (pt3, pt4);
    double distance_2D_ans = std::sqrt(18);

    if (std::fabs(distance_2D - distance_2D_ans) > SMALLNUMBER)
    {
        std::cout << "test_function_8: distance 2D failed" << std::endl;
        std::cout << "distance_2D: " << distance_2D << std::endl;
        std::cout << "distance_2D_ans: " << distance_2D_ans << std::endl;
        return false;
    }

    return true;
}

// test calculate the distance between point and line
bool test_function_9 ()
{
    Pt3D pt1(1,2,3);
    Pt3D pt2(4,5,6);
    Line3D line;
    line.pt = pt1;
    line.unit_vector = myMATH::createUnitVector (pt1, pt2);
    Pt3D pt3(2,3,4);
    double distance_3D = myMATH::dist (pt3, line);
    double distance_3D_ans = 0;

    if (std::fabs(distance_3D - distance_3D_ans) > SQRTSMALLNUMBER)
    {
        std::cout << "test_function_9: distance 3D failed" << std::endl;
        std::cout << "distance_3D: " << distance_3D << std::endl;
        std::cout << "distance_3D_ans: " << distance_3D_ans << std::endl;
        line.pt.print();
        line.unit_vector.print();
        return false;
    }

    Pt2D pt4(1,2);
    Pt2D pt5(4,5);
    Line2D line_2D;
    line_2D.pt = pt4;
    line_2D.unit_vector = myMATH::createUnitVector (pt4, pt5);
    Pt2D pt6(10,13);
    double distance_2D = myMATH::dist (pt6, line_2D);
    double distance_2D_ans = std::fabs(pt6[0] - pt6[1] + 1) / std::sqrt(2);

    if (std::fabs(distance_2D - distance_2D_ans) > SMALLNUMBER)
    {
        std::cout << "test_function_9: distance 2D failed" << std::endl;
        std::cout << "distance_2D: " << distance_2D << std::endl;
        std::cout << "distance_2D_ans: " << distance_2D_ans << std::endl;
        return false;
    }

    return true;
}

// test create identity matrix
bool test_function_10 ()
{
    Matrix<double> mtx = myMATH::eye<double> (3);
    Matrix<double> mtx_ans({{1,0,0},{0,1,0},{0,0,1}});

    if (mtx != mtx_ans)
    {
        std::cout << "test_function_10: createIdentityMatrix failed" << std::endl;
        std::cout << "mtx: " << std::endl;
        mtx.print();
        std::cout << "mtx_ans: " << std::endl;
        mtx_ans.print();
        return false;
    }

    return true;
}

// test piecewise product
bool test_function_11 ()
{
    Matrix<double> mtx1({{1,2,3},{4,5,6},{7,8,9}});
    Matrix<double> mtx2({{9,8,7},{6,5,4},{3,2,1}});
    Matrix<double> mtx = myMATH::piecewiseProduct (mtx1, mtx2);
    Matrix<double> mtx_ans({{9,16,21},{24,25,24},{21,16,9}});

    if (mtx != mtx_ans)
    {
        std::cout << "test_function_11: piecewiseProduct failed" << std::endl;
        std::cout << "mtx: " << std::endl;
        mtx.print();
        std::cout << "mtx_ans: " << std::endl;
        mtx_ans.print();
        return false;
    }

    return true;
}

// test matrix inverse
bool test_function_12 ()
{
    Matrix<double> mtx({{1,2,3},{4,5,6},{7,8,10}});
    
    Matrix<double> mtx_inv_gauss = myMATH::inverse (mtx, "gauss");
    Matrix<double> mtx_inv_det = myMATH::inverse (mtx, "det");
    Matrix<double> mtx_inv_ans(
       {{-0.66666667, -1.33333333,  1.},
        {-0.66666667,  3.66666667, -2.},
        { 1.        , -2.        ,  1.}});

    if (mtx_inv_gauss != mtx_inv_ans)
    {
        std::cout << "test_function_12: matrix inverse (gauss) failed" << std::endl;
        std::cout << "mtx_inv: " << std::endl;
        mtx_inv_gauss.print();
        std::cout << "mtx_inv_ans: " << std::endl;
        mtx_inv_ans.print();
        return false;
    }
    if (mtx_inv_det != mtx_inv_ans)
    {
        std::cout << "test_function_12: matrix inverse (det) failed" << std::endl;
        std::cout << "mtx_inv: " << std::endl;
        mtx_inv_det.print();
        std::cout << "mtx_inv_ans: " << std::endl;
        mtx_inv_ans.print();
        return false;
    }

    Matrix<int> mtx_2({{1,2,3},{4,5,6},{7,8,9}});
    try
    {
        Matrix<int> mtx_inv = myMATH::inverse (mtx_2, "gauss");
    }
    catch (const std::exception& e)
    {
        std::cout << "test_function_12: matrix inverse (gauss) failed" << std::endl;
        std::cout << e.what() << std::endl;
        return false;
    }
    
    // test 2x2 matrix
    Matrix<double> mtx_3({{1,2},{3,4}});
    mtx_inv_gauss = myMATH::inverse (mtx_3, "gauss");
    mtx_inv_det = myMATH::inverse (mtx_3, "det");
    mtx_inv_ans = Matrix<double>({{-2,1},{1.5,-0.5}});

    if (mtx_inv_gauss != mtx_inv_ans)
    {
        std::cout << "test_function_12 (line " << __LINE__ 
                  << "): matrix inverse (gauss) failed" << std::endl;
        std::cout << "mtx_inv: " << std::endl;
        mtx_inv_gauss.print();
        std::cout << "mtx_inv_ans: " << std::endl;
        mtx_inv_ans.print();
        return false;
    }
    if (mtx_inv_det != mtx_inv_ans)
    {
        std::cout << "test_function_12 (line " << __LINE__ 
                  << "): matrix inverse (det) failed" << std::endl;
        std::cout << "mtx_inv: " << std::endl;
        mtx_inv_det.print();
        std::cout << "mtx_inv_ans: " << std::endl;
        mtx_inv_ans.print();
        return false;
    }

    return true;
}

// test matrix trace
bool test_function_13 ()
{
    Matrix<double> mtx({{1,2,3},{4,5,6},{7,8,9}});
    double trace = myMATH::trace (mtx);
    double trace_ans = 15;

    if (std::fabs(trace - trace_ans) > SMALLNUMBER)
    {
        std::cout << "test_function_13: matrix trace failed" << std::endl;
        std::cout << "trace: " << trace << std::endl;
        std::cout << "trace_ans: " << trace_ans << std::endl;
        return false;
    }

    return true;
}

// test triangulation 
bool test_function_14 ()
{
    // CASE 1
    std::vector<Line3D> line_of_sight_list;
    Line3D line1;
    line1.pt = Pt3D(0,0,0);
    line1.unit_vector = Pt3D(1,0,0);
    line_of_sight_list.push_back(line1);
    Line3D line2;
    line2.pt = Pt3D(0,0,0);
    line2.unit_vector = Pt3D(0,1,0);
    line_of_sight_list.push_back(line2);
    Line3D line3;
    line3.pt = Pt3D(0,0,0);
    line3.unit_vector = Pt3D(0,0,1);
    line_of_sight_list.push_back(line3);

    Pt3D pt_world(1,2,3);
    double error;

    myMATH::triangulation (pt_world, error, line_of_sight_list);

    if (pt_world != Pt3D(0,0,0) || error > SMALLNUMBER)
    {
        std::cout << "test_function_14: triangulation failed at line " << __LINE__ << std::endl;
        std::cout << "pt_world: " << std::endl;
        pt_world.print();
        std::cout << "pt_world_ans: " << std::endl;
        Pt3D(0,0,0).print();
        return false;
    }

    // CASE 2
    line_of_sight_list.clear();
    line1.pt = Pt3D(0,0,0);
    line1.unit_vector = Pt3D(1,1,1) / std::sqrt(3);
    line_of_sight_list.push_back(line1);

    line2.pt = Pt3D(1,0,0);
    line2.unit_vector = Pt3D(1,-1,-1) / std::sqrt(3);
    line_of_sight_list.push_back(line2);

    myMATH::triangulation (pt_world, error, line_of_sight_list);

    if (pt_world != Pt3D(0.5,0.5,0.5) || error > SMALLNUMBER)
    {
        std::cout << "test_function_14: triangulation failed at line " << __LINE__ << std::endl;
        std::cout << "pt_world: " << std::endl;
        pt_world.print();
        return false;
    }

    return true;
}

// test isLocalMax
bool test_function_15 ()
{
    Matrix<double> img({{1,0,0,0,0},
                        {0,1,2,1,0},
                        {0,2,3,2,0},
                        {0,1,2,1,0},
                        {0,0,0,0,0}});
    bool is_localmax = myMATH::isLocalMax (img, 2, 2);
    bool is_localmax_ans = true;
    if (is_localmax != is_localmax_ans)
    {
        std::cout << "test_function_15: isLocalMax failed at line " << __LINE__ << std::endl;
        std::cout << "is_localmax: " << is_localmax << std::endl;
        std::cout << "is_localmax_ans: " << is_localmax_ans << std::endl;
        return false;
    }

    is_localmax = myMATH::isLocalMax (img, 1, 1);
    is_localmax_ans = false;
    if (is_localmax != is_localmax_ans)
    {
        std::cout << "test_function_15: isLocalMax failed at line" << __LINE__ << std::endl;
        std::cout << "is_localmax: " << is_localmax << std::endl;
        std::cout << "is_localmax_ans: " << is_localmax_ans << std::endl;
        return false;
    }

    is_localmax = myMATH::isLocalMax (img, 0, 0);
    is_localmax_ans = true;
    if (is_localmax != is_localmax_ans)
    {
        std::cout << "test_function_15: isLocalMax failed at line" << __LINE__ << std::endl;
        std::cout << "is_localmax: " << is_localmax << std::endl;
        std::cout << "is_localmax_ans: " << is_localmax_ans << std::endl;
        return false;
    }

    return true;
}

// test find cross point of two 2d lines
bool test_function_16 ()
{
    Pt2D pt1(0,1);
    Pt2D pt2(3,4);
    Line2D line1;
    line1.pt = pt1;
    line1.unit_vector = myMATH::createUnitVector (pt1, pt2);
    
    Pt2D pt3(1,2);
    Pt2D pt4(2,1);
    Line2D line2;
    line2.pt = pt3;
    line2.unit_vector = myMATH::createUnitVector (pt3, pt4);

    Pt2D cross_point = myMATH::crossPoint (line1, line2);

    if (cross_point != Pt2D(1,2))
    {
        std::cout << "test_function_16: crossPoint failed" << std::endl;
        std::cout << "cross_point: " << std::endl;
        cross_point.print();
        std::cout << "cross_point_ans: " << std::endl;
        Pt2D(1,2).print();
        return false;
    }



    return true;
}

// test polyfit
bool test_function_17 ()
{
    // load x_list
    std::vector<std::vector<double>> x_list;
    std::ifstream file("../test/inputs/test_myMATH/x_list.txt");
    std::string line;
    std::cout << "load x:" << std::endl;
    while (std::getline(file, line))
    {
        std::stringstream is(line);
        std::vector<double> val_list;
        std::string val;
        while (std::getline(is, val, ','))
        {
            val_list.push_back(std::stof(val));
        }
        x_list.push_back(val_list);
    }
    file.close();

    // load y_list
    std::vector<std::vector<double>> y_list;
    file.open("../test/inputs/test_myMATH/y_list.txt");
    while (std::getline(file, line))
    {
        std::stringstream is(line);
        std::vector<double> val_list;
        std::string val;
        while (std::getline(is, val, ','))
        {
            val_list.push_back(std::stof(val));
        }
        y_list.push_back(val_list);
    }
    file.close();

    // load order_list
    std::vector<int> order_list;
    file.open("../test/inputs/test_myMATH/order_list.txt");
    while (std::getline(file, line))
    {
        order_list.push_back(std::stoi(line));
    }
    file.close();

    // load sol: coeff_list
    std::vector<std::vector<double>> coeff_list_sol;
    file.open("../test/solutions/test_myMATH/coeff_list.txt");
    while (std::getline(file, line))
    {
        std::stringstream is(line);
        std::vector<double> val_list;
        std::string val;
        while (std::getline(is, val, ','))
        {
            val_list.push_back(std::stof(val));
        }
        coeff_list_sol.push_back(val_list);
    }
    file.close();

    // calculate coeff_list using polyfit
    std::vector<std::vector<double>> coeff_list;
    for (int i = 0; i < x_list.size(); i ++)
    {
        std::vector<double> coeff;
        myMATH::polyfit (coeff, x_list[i], y_list[i], order_list[i]);
        coeff_list.push_back(coeff);
    }

    for (int i = 0; i < coeff_list_sol.size(); i ++)
    {
        for (int j = 0; j < coeff_list_sol[i].size(); j ++)
        {
            if (std::fabs(coeff_list_sol[i][j]-coeff_list[i][j]) > 1e-5)
            {
                std::cout << "test_function_17: polyfit failed" << std::endl;
                std::cout << "coeff: " << std::endl;
                for (int k = 0; k < coeff_list[i].size(); k ++)
                {
                    std::cout << std::setprecision(SAVEPRECISION) << coeff_list[i][k] << " ";
                }
                std::cout << std::endl;

                std::cout << "coeff_ans: " << std::endl;
                for (int k = 0; k < coeff_list_sol[i].size(); k ++)
                {
                    std::cout << std::setprecision(SAVEPRECISION) << coeff_list_sol[i][k] << " ";
                }
                std::cout << std::endl;
                
                return false;
            }
        }
    }

    return true;
}

int main()
{
    fs::create_directories("../test/results/test_myMATH/");

    IS_TRUE(test_function_1());
    IS_TRUE(test_function_2());
    IS_TRUE(test_function_3());
    IS_TRUE(test_function_4());
    IS_TRUE(test_function_5());
    IS_TRUE(test_function_6());
    IS_TRUE(test_function_7());
    IS_TRUE(test_function_8());
    IS_TRUE(test_function_9());
    IS_TRUE(test_function_10());
    IS_TRUE(test_function_11());
    IS_TRUE(test_function_12());
    IS_TRUE(test_function_13());
    IS_TRUE(test_function_14());
    IS_TRUE(test_function_15());
    IS_TRUE(test_function_16());
    IS_TRUE(test_function_17());

    return 0;
}