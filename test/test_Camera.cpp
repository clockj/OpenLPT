#include "test.h"
#include "myMATH.h"
#include "Matrix.h"
#include "Camera.h"

// pinhole camera

// test constructor
bool test_function_1 ()
{
    Camera c;
    for (int i = 0; i < 5; i ++)
    {
        c.loadParameters("../test/inputs/test_Camera/cam"+std::to_string(i+1)+".txt");
        c.saveParameters("../test/results/test_Camera/cam"+std::to_string(i+1)+".txt");
    }
    
    // compare the results with the inputs
    for (int i = 0; i < 5; i ++)
    {
        std::ifstream file1("../test/results/test_Camera/cam"+std::to_string(i+1)+".txt");
        std::ifstream file2("../test/solutions/test_Camera/cam"+std::to_string(i+1)+".txt");

        std::string line1, line2;
        std::stringstream ss1, ss2;
        while (std::getline(file1, line1) && std::getline(file2, line2))
        {
            if (line1 != line2)
            {
                std::stringstream ss1, ss2;
                ss1 << line1;
                ss2 << line2;
                double a, b;
                while (ss1 >> a && ss2 >> b)
                {
                    // std::cout << a << "," << b << std::endl;
                    if (std::fabs(a - b) > 1e-5)
                    {
                        std::cout << "test_function_1: failed at cam" << i+1 << std::endl;
                        return false;
                    }

                    if (ss1.peek() == ',')
                    {
                        ss1.ignore();
                    }
                    if (ss2.peek() == ',')
                    {
                        ss2.ignore();
                    }
                }
            }
        }
    }
    
    return true;
}

// test projection
bool test_function_2 ()
{
    Pt3D pt_world(1, 2, 3);
    Camera c;

    // load solutions
    std::ifstream file("../test/solutions/test_Camera/test_function_2.txt");
    std::string line;
    std::stringstream is;
    while (std::getline(file, line))
    {
        is << line << '\t';
    }
    file.close();

    // test worldToUndistImg
    for (int i = 0; i < 5; i ++)
    {
        c.loadParameters("../test/inputs/test_Camera/cam"+std::to_string(i+1)+".txt");
        Pt2D pt_undist = c.worldToUndistImg(pt_world);
        Pt2D pt_solution(is);
        
        if (pt_undist != pt_solution)
        {
            std::cout << "test_function_2 (line " << __LINE__ << "): failed at cam" << i+1 << std::endl;
            pt_undist.print(8);
            pt_solution.print(8);
            return false;
        }
    }

    // test distortion
    for (int i = 0; i < 5; i ++)
    {
        c.loadParameters("../test/inputs/test_Camera/cam"+std::to_string(i+1)+".txt");
        Pt2D pt_undist = c.worldToUndistImg(pt_world);
        Pt2D pt_dist = c.distort(pt_undist);
        Pt2D pt_solution(is);
        
        // tolerance set to be 1e-3
        Pt2D residue = pt_dist - pt_solution;
        if (std::fabs(residue[0]) > 1e-3 || std::fabs(residue[1]) > 1e-3)
        {
            std::cout << "test_function_2 (line " << __LINE__ << "): failed at cam" << i+1 << std::endl;
            pt_dist.print(8);
            pt_solution.print(8);
            return false;
        }
    }

    // test project
    for (int i = 0; i < 5; i ++)
    {
        c.loadParameters("../test/inputs/test_Camera/cam"+std::to_string(i+1)+".txt");
        Pt2D pt_dist = c.project(pt_world);
        Pt2D pt_solution(is);
        
        // tolerance set to be 1e-3
        Pt2D residue = pt_dist - pt_solution;
        if (std::fabs(residue[0]) > 1e-3 || std::fabs(residue[1]) > 1e-3)
        {
            std::cout << "test_function_2 (line " << __LINE__ << "): failed at cam" << i+1 << std::endl;
            pt_dist.print(8);
            pt_solution.print(8);
            return false;
        }
    }

    return true;
}

// test line of sight 
bool test_function_3 ()
{
    Camera c;
    Pt2D pt_dist(461, 441);

    // load solutions
    std::ifstream file("../test/solutions/test_Camera/test_function_3.txt");
    std::string line;
    std::stringstream is;
    while (std::getline(file, line))
    {
        is << line << '\t';
    }
    file.close();

    // test undistort
    for (int i = 0; i < 5; i ++)
    {
        c.loadParameters("../test/inputs/test_Camera/cam"+std::to_string(i+1)+".txt");
        Pt2D pt_undist = c.undistort(pt_dist);
        Pt2D pt_sol(is);
        
        // tolerance set to be 1e-5
        Pt2D pt_residue = pt_undist - pt_sol;
        if (std::fabs(pt_residue[0]) > 1e-5 || std::fabs(pt_residue[1]) > 1e-5)
        {
            std::cout << "test_function_3 (line " << __LINE__ << "): failed at cam" << i+1 << std::endl;
            pt_undist.print(8);
            pt_sol.print(8);
            return false;
        }
    }

    // test pinholeLine
    for (int i = 0; i < 5; i ++)
    {
        c.loadParameters("../test/inputs/test_Camera/cam"+std::to_string(i+1)+".txt");
        Line3D line = c.pinholeLine(c.undistort(pt_dist));
        
        Pt3D pt_ref(is);
        Pt3D unit_vec(is);
        
        // tolerance set to be 1e-5
        Pt3D pt_residue = line.pt - pt_ref;
        Pt3D vec_residue = line.unit_vector - unit_vec;
        for (int j = 0; j < 3; j ++)
        {
            if (std::fabs(pt_residue[j]) > 1e-5 || std::fabs(vec_residue[j]) > 1e-5)
            {
                std::cout << "test_function_3 (line " << __LINE__ << "): failed at cam" << i+1 << std::endl;
                line.pt.print(8);
                line.unit_vector.print(8);

                pt_ref.print(8);
                unit_vec.print(8);
                return false;
            }
        }
    }

    // test lineOfSight
    for (int i = 0; i < 5; i ++)
    {
        c.loadParameters("../test/inputs/test_Camera/cam"+std::to_string(i+1)+".txt");
        Line3D line = c.lineOfSight(pt_dist);
        
        Pt3D pt_ref(is);
        Pt3D unit_vec(is);
        
        // tolerance set to be 1e-5
        Pt3D pt_residue = line.pt - pt_ref;
        Pt3D vec_residue = line.unit_vector - unit_vec;
        for (int j = 0; j < 3; j ++)
        {
            if (std::fabs(pt_residue[j]) > 1e-5 || std::fabs(vec_residue[j]) > 1e-5)
            {
                std::cout << "test_function_3 (line " << __LINE__ << "): failed at cam" << i+1 << std::endl;
                line.pt.print(8);
                line.unit_vector.print(8);

                pt_ref.print(8);
                unit_vec.print(8);
                return false;
            }
        }
    }

    return true;
}


// Polynomial model

// test constructor
bool test_function_4 ()
{
    Camera c;
    for (int i = 0; i < 2; i ++)
    {
        c.loadParameters("../test/inputs/test_Camera/cam"+std::to_string(i+1)+"_poly"+".txt");
        c.saveParameters("../test/results/test_Camera/cam"+std::to_string(i+1)+"_poly"+".txt");

        c._poly_param.du_coeffs.write("../test/results/test_Camera/cam"+std::to_string(i+1)+"_poly_du"+".txt");
        c._poly_param.dv_coeffs.write("../test/results/test_Camera/cam"+std::to_string(i+1)+"_poly_dv"+".txt");
    }
    
    // compare the results with the inputs
    for (int i = 0; i < 2; i ++)
    {
        std::ifstream file1("../test/results/test_Camera/cam"+std::to_string(i+1)+"_poly"+".txt");
        std::ifstream file2("../test/solutions/test_Camera/cam"+std::to_string(i+1)+"_poly"+".txt");

        std::string line1, line2;
        std::stringstream ss1, ss2;
        while (std::getline(file1, line1) && std::getline(file2, line2))
        {
            if (line1 != line2)
            {
                std::cout << "test_function_4: failed at cam" << i+1 << std::endl;
                std::cout << line1 << std::endl;
                std::cout << line2 << std::endl;
                return false;
            }
        }
    }
    
    return true;
}

// test projection
bool test_function_5 ()
{
    Pt3D pt_world(1, 0, 3);
    Camera c;

    // load solutions
    std::ifstream file("../test/solutions/test_Camera/test_function_5.txt");
    std::string line;
    std::stringstream is;
    while (std::getline(file, line))
    {
        is << line << '\t';
    }
    file.close();

    // test polyProject
    for (int i = 0; i < 2; i ++)
    {
        c.loadParameters("../test/inputs/test_Camera/cam"+std::to_string(i+1)+"_poly"+".txt");
        Pt2D pt_dist = c.polyProject(pt_world);
        Pt2D pt_solution(is);
        
        if (pt_dist != pt_solution)
        {
            std::cout << "test_function_5 (line " << __LINE__ << "): failed at cam" << i+1 << std::endl;
            pt_dist.print(8);
            pt_solution.print(8);
            return false;
        }
    }

    // test project
    for (int i = 0; i < 2; i ++)
    {
        c.loadParameters("../test/inputs/test_Camera/cam"+std::to_string(i+1)+"_poly"+".txt");
        Pt2D pt_dist = c.project(pt_world);
        Pt2D pt_solution(is);
        
        // tolerance set to be 1e-3
        if (pt_dist != pt_solution)
        {
            std::cout << "test_function_5 (line " << __LINE__ << "): failed at cam" << i+1 << std::endl;
            pt_dist.print(8);
            pt_solution.print(8);
            return false;
        }
    }

    return true;
}

// test line of sight
bool test_function_6 ()
{
    Pt3D pt_world(1, 0, 3);
    Camera c;

    for (int i = 0; i < 4; i ++)
    {
        c.loadParameters("../test/inputs/test_Camera/cam"+std::to_string(i+1)+"_poly"+".txt");
        c.saveParameters("../test/results/test_Camera/cam"+std::to_string(i+1)+"_poly"+".txt");

        Pt2D pt_dist = c.polyProject(pt_world);
        Line3D line = c.polyLineOfSight(pt_dist);

        // check whether the orginal point is on the line
        Pt3D vec = pt_world - line.pt;
        double cos = myMATH::dot(vec, line.unit_vector) / (vec.norm() * line.unit_vector.norm());
        std::cout << "cos = " << cos << std::endl;
        if (std::fabs(cos) < 0.98 || std::isnan(cos))
        {
            std::cout << "test_function_6 (line " << __LINE__ << "): failed at cam" << i+1 << ": cos = " << cos << std::endl;

            pt_world.print(8);
            line.pt.print(8);
            line.unit_vector.print(8);
            
            return false;
        }
    }
    
    for (int i = 0; i < 4; i ++)
    {
        c.loadParameters("../test/inputs/test_Camera/cam"+std::to_string(i+1)+"_poly"+".txt");

        Pt2D pt_dist = c.project(pt_world);
        Line3D line = c.lineOfSight(pt_dist);

        // check whether the orginal point is on the line
        Pt3D vec = pt_world - line.pt;
        double cos = myMATH::dot(vec, line.unit_vector) / (vec.norm() * line.unit_vector.norm());
        std::cout << "cos = " << cos << std::endl;
        if (std::fabs(cos) < 0.98 || std::isnan(cos))
        {
            std::cout << "test_function_6 (line " << __LINE__ << "): failed at cam" << i+1 << ": cos = " << cos << std::endl;

            pt_world.print(8);
            line.pt.print(8);
            line.unit_vector.print(8);
            
            return false;
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

    return 0;
}

