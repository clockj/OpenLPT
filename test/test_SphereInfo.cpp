#include "test.h"
#include "SphereInfo.h"
#include "Matrix.h"
#include "Camera.h"

// test constructor
bool test_function_1 ()
{
    Circle circle(Pt2D(3, 4));
    if (circle._pt_center != Pt2D(3, 4))
    {
        return false;
    }

    Sphere sphere(Pt3D(8, 9, 10));
    if (sphere._pt_center != Pt3D(8, 9, 10))
    {
        return false;
    }

    return true;
}

// test Sphere functions
bool test_function_2 ()
{
    Sphere sphere;
    sphere._pt_center = Pt3D(1, 2, 3);
    std::vector<int> camid_list;
    std::vector<Circle> circle_list;

    // test addCircle
    Circle circle_1 = Circle(Pt2D(4, 5));
    Circle circle_2 = Circle(Pt2D(6, 7));
    Circle circle_3 = Circle(Pt2D(8, 9));
    std::vector<Circle> circle_list_ans = {circle_1, circle_2, circle_3};
    std::vector<int> camid_list_ans = {0, 10, 2};

    sphere.addCircle(circle_1, 0);
    sphere.addCircle(circle_2, 10);
    sphere.addCircle(circle_3, 2);
    camid_list = sphere._camid_list;
    circle_list = sphere._circle_list;
    if (camid_list.size() != 3 || circle_list.size() != 3)
    {
        std::cout << "test_function_2: addCircle failed (line " << __LINE__ << ")" << std::endl;
        std::cout << "camid_list.size() = " << camid_list.size() << std::endl;
        std::cout << "circle_list.size() = " << circle_list.size() << std::endl;
        return false;
    }
    for (int i = 0; i < 3; i++)
    {
        if (camid_list[i] != camid_list_ans[i] || circle_list[i]._pt_center != circle_list_ans[i]._pt_center)
        {
            std::cout << "test_function_2: addCircle failed (line " << __LINE__ << ")" << std::endl;
            std::cout << "camid_list[i] = " << camid_list[i] << std::endl;
            std::cout << "camid_list_ans[i] = " << camid_list_ans[i] << std::endl;
            return false;
        }
    }
    ////////////////////////

    // test removeCircle
    sphere.removeCircle(100);
    camid_list = sphere._camid_list;
    circle_list = sphere._circle_list;
    if (camid_list.size() != 3 || circle_list.size() != 3)
    {
        std::cout << "test_function_2: removeCircle failed (line " << __LINE__ << ")" << std::endl;
        std::cout << "camid_list.size() = " << camid_list.size() << std::endl;
        std::cout << "circle_list.size() = " << circle_list.size() << std::endl;
        return false;
    }
    for (int i = 0; i < 3; i++)
    {
        if (camid_list[i] != camid_list_ans[i] || circle_list[i]._pt_center != circle_list_ans[i]._pt_center)
        {
            std::cout << "test_function_2: removeCircle failed (line " << __LINE__ << ")" << std::endl;
            std::cout << "camid_list[i] = " << camid_list[i] << std::endl;
            std::cout << "camid_list_ans[i] = " << camid_list_ans[i] << std::endl;
            return false;
        }
    }

    sphere.removeCircle(10);
    std::vector<Circle> circle_list_ans_2 = {circle_1, circle_3};
    std::vector<int> camid_list_ans_2 = {0, 2};
    camid_list = sphere._camid_list;
    circle_list = sphere._circle_list;
    if (camid_list.size() != 2 || circle_list.size() != 2)
    {
        std::cout << "test_function_2: removeCircle failed (line " << __LINE__ << ")" << std::endl;
        std::cout << "camid_list.size() = " << camid_list.size() << std::endl;
        std::cout << "circle_list.size() = " << circle_list.size() << std::endl;
        return false;
    }
    for (int i = 0; i < 2; i++)
    {
        if (camid_list[i] != camid_list_ans_2[i] || circle_list[i]._pt_center != circle_list_ans_2[i]._pt_center)
        {
            std::cout << "test_function_2: removeCircle failed (line " << __LINE__ << ")" << std::endl;
            std::cout << "camid_list[i] = " << camid_list[i] << std::endl;
            std::cout << "camid_list_ans[i] = " << camid_list_ans[i] << std::endl;
            return false;
        }
    }

    sphere.removeCircle({0, 2});
    camid_list = sphere._camid_list;
    circle_list = sphere._circle_list;
    if (camid_list.size() != 0 || circle_list.size() != 0)
    {
        std::cout << "test_function_2: removeCircle failed (line " << __LINE__ << ")" << std::endl;
        std::cout << "camid_list.size() = " << camid_list.size() << std::endl;
        std::cout << "circle_list.size() = " << circle_list.size() << std::endl;
        return false;
    }
    ////////////////////////

    // test clearCircle
    sphere.addCircle(circle_1, 0);
    sphere.addCircle(circle_2, 10);
    sphere.addCircle(circle_3, 2);
    sphere.clearCircle();
    camid_list = sphere._camid_list;
    circle_list = sphere._circle_list;
    if (camid_list.size() != 0 || circle_list.size() != 0)
    {
        std::cout << "test_function_2: clearCircle failed (line " << __LINE__ << ")" << std::endl;
        std::cout << "camid_list.size() = " << camid_list.size() << std::endl;
        std::cout << "circle_list.size() = " << circle_list.size() << std::endl;
        return false;
    }
    ////////////////////////

    // test updateCircle
    sphere.addCircle(circle_1, 0);
    sphere.addCircle(circle_2, 10);
    sphere.addCircle(circle_3, 2);
    Circle circle_4 = Circle(Pt2D(11, 12));
    sphere.updateCircle(circle_4, 10);
    circle_list = sphere._circle_list;
    if (circle_list.size() != 3 || circle_list[1]._pt_center != circle_4._pt_center)
    {
        std::cout << "test_function_2: updateCircle failed (line " << __LINE__ << ")" << std::endl;
        std::cout << "circle_list.size() = " << circle_list.size() << std::endl;
        return false;
    }

    sphere.updateCircle({circle_4, circle_4}, {0, 2});
    circle_list = sphere._circle_list;
    if (circle_list.size() != 2 || 
        circle_list[0]._pt_center != circle_4._pt_center || 
        circle_list[1]._pt_center != circle_4._pt_center)
    {
        std::cout << "test_function_2: updateCircle failed (line " << __LINE__ << ")" << std::endl;
        std::cout << "circle_list.size() = " << circle_list.size() << std::endl;
        return false;
    }
    ////////////////////////

    // test projectObject2D
    std::vector<Camera> cam_list_all;
    std::vector<int> camid_list_2 = {0, 1, 2, 3};
    for (int i = 0; i < 5; i ++)
    {
        Camera cam ("../test/inputs/test_ObjectInfo/cam"+std::to_string(i+1)+".txt");
        cam_list_all.push_back(cam);
    }
    sphere.projectObject2D(camid_list_2, cam_list_all);
    circle_list = sphere._circle_list;
    camid_list = sphere._camid_list;
    if (circle_list.size() != 4 || camid_list.size() != 4)
    {
        std::cout << "test_function_2: projectObject2D failed (line " << __LINE__ << ")" << std::endl;
        std::cout << "circle_list.size() = " << circle_list.size() << std::endl;
        std::cout << "camid_list.size() = " << camid_list.size() << std::endl;
        return false;
    }
    for (int i = 0; i < 4; i++)
    {
        if (circle_list[i]._pt_center != cam_list_all[i].project(sphere._pt_center) || camid_list[i] != camid_list_2[i])
        {
            std::cout << "test_function_2: projectObject2D failed (line " << __LINE__ << ")" << std::endl;
            return false;
        }
    }
    ////////////////////////

    return true;
}

int main ()
{
    fs::create_directory("../test/results/test_SphereInfo");

    IS_TRUE(test_function_1());
    IS_TRUE(test_function_2());

    return 0;
}