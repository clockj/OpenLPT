#include "test.h"
#include "ObjectInfo.h"
#include "Camera.h"

// test constructor
bool test_function_1 ()
{
    Object2D obj2d(Pt2D(1, 2));
    if (obj2d._pt_center != Pt2D(1, 2))
    {
        return false;
    }

    Tracer2D tr2d(Pt2D(3, 4));
    if (tr2d._pt_center != Pt2D(3, 4))
    {
        return false;
    }

    Object3D obj3d(Pt3D(5, 6, 7));
    if (obj3d._pt_center != Pt3D(5, 6, 7))
    {
        return false;
    }

    Tracer3D tr3d(Pt3D(8, 9, 10));
    if (tr3d._pt_center != Pt3D(8, 9, 10))
    {
        return false;
    }


    return true;
}

// test Tracer3D functions
bool test_function_2 ()
{
    Tracer3D tr3d;
    tr3d._pt_center = Pt3D(1, 2, 3);
    std::vector<int> camid_list;
    std::vector<Tracer2D> tracer2d_list;

    // test addTracer2D
    Tracer2D tr2d_1 = Tracer2D(Pt2D(4, 5));
    Tracer2D tr2d_2 = Tracer2D(Pt2D(6, 7));
    Tracer2D tr2d_3 = Tracer2D(Pt2D(8, 9));
    std::vector<Tracer2D> tr2d_list_ans = {tr2d_1, tr2d_2, tr2d_3};
    std::vector<int> camid_list_ans = {0, 10, 2};

    tr3d.addTracer2D(tr2d_1, 0);
    tr3d.addTracer2D(tr2d_2, 10);
    tr3d.addTracer2D(tr2d_3, 2);
    camid_list = tr3d._camid_list;
    tracer2d_list = tr3d._tr2d_list;
    if (camid_list.size() != 3 || tracer2d_list.size() != 3)
    {
        std::cout << "test_function_2: addTracer2D failed (line " << __LINE__ << ")" << std::endl;
        std::cout << "camid_list.size() = " << camid_list.size() << std::endl;
        std::cout << "tracer2d_list.size() = " << tracer2d_list.size() << std::endl;
        return false;
    }
    for (int i = 0; i < 3; i++)
    {
        if (camid_list[i] != camid_list_ans[i] || tracer2d_list[i]._pt_center != tr2d_list_ans[i]._pt_center)
        {
            std::cout << "test_function_2: addTracer2D failed (line " << __LINE__ << ")" << std::endl;
            std::cout << "camid_list[i] = " << camid_list[i] << std::endl;
            std::cout << "camid_list_ans[i] = " << camid_list_ans[i] << std::endl;
            return false;
        }
    }
    ////////////////////////

    // test removeTracer2D
    tr3d.removeTracer2D(100);
    camid_list = tr3d._camid_list;
    tracer2d_list = tr3d._tr2d_list;
    if (camid_list.size() != 3 || tracer2d_list.size() != 3)
    {
        std::cout << "test_function_2: removeTracer2D failed (line " << __LINE__ << ")" << std::endl;
        std::cout << "camid_list.size() = " << camid_list.size() << std::endl;
        std::cout << "tracer2d_list.size() = " << tracer2d_list.size() << std::endl;
        return false;
    }
    for (int i = 0; i < 3; i++)
    {
        if (camid_list[i] != camid_list_ans[i] || tracer2d_list[i]._pt_center != tr2d_list_ans[i]._pt_center)
        {
            std::cout << "test_function_2: removeTracer2D failed (line " << __LINE__ << ")" << std::endl;
            std::cout << "camid_list[i] = " << camid_list[i] << std::endl;
            std::cout << "camid_list_ans[i] = " << camid_list_ans[i] << std::endl;
            return false;
        }
    }

    tr3d.removeTracer2D(10);
    std::vector<Tracer2D> tr2d_list_ans_2 = {tr2d_1, tr2d_3};
    std::vector<int> camid_list_ans_2 = {0, 2};
    camid_list = tr3d._camid_list;
    tracer2d_list = tr3d._tr2d_list;
    if (camid_list.size() != 2 || tracer2d_list.size() != 2)
    {
        std::cout << "test_function_2: removeTracer2D failed (line " << __LINE__ << ")" << std::endl;
        std::cout << "camid_list.size() = " << camid_list.size() << std::endl;
        std::cout << "tracer2d_list.size() = " << tracer2d_list.size() << std::endl;
        return false;
    }
    for (int i = 0; i < 2; i++)
    {
        if (camid_list[i] != camid_list_ans_2[i] || tracer2d_list[i]._pt_center != tr2d_list_ans_2[i]._pt_center)
        {
            std::cout << "test_function_2: removeTracer2D failed (line " << __LINE__ << ")" << std::endl;
            std::cout << "camid_list[i] = " << camid_list[i] << std::endl;
            std::cout << "camid_list_ans[i] = " << camid_list_ans[i] << std::endl;
            return false;
        }
    }

    tr3d.removeTracer2D({0, 2});
    camid_list = tr3d._camid_list;
    tracer2d_list = tr3d._tr2d_list;
    if (camid_list.size() != 0 || tracer2d_list.size() != 0)
    {
        std::cout << "test_function_2: removeTracer2D failed (line " << __LINE__ << ")" << std::endl;
        std::cout << "camid_list.size() = " << camid_list.size() << std::endl;
        std::cout << "tracer2d_list.size() = " << tracer2d_list.size() << std::endl;
        return false;
    }
    ////////////////////////

    // test clearTracer2D
    tr3d.addTracer2D(tr2d_1, 0);
    tr3d.addTracer2D(tr2d_2, 10);
    tr3d.addTracer2D(tr2d_3, 2);
    tr3d.clearTracer2D();
    camid_list = tr3d._camid_list;
    tracer2d_list = tr3d._tr2d_list;
    if (camid_list.size() != 0 || tracer2d_list.size() != 0)
    {
        std::cout << "test_function_2: clearTracer2D failed (line " << __LINE__ << ")" << std::endl;
        std::cout << "camid_list.size() = " << camid_list.size() << std::endl;
        std::cout << "tracer2d_list.size() = " << tracer2d_list.size() << std::endl;
        return false;
    }
    ////////////////////////

    // test updateTracer2D
    tr3d.addTracer2D(tr2d_1, 0);
    tr3d.addTracer2D(tr2d_2, 10);
    tr3d.addTracer2D(tr2d_3, 2);
    Tracer2D tr2d_4 = Tracer2D(Pt2D(11, 12));
    tr3d.updateTracer2D(tr2d_4, 10);
    tracer2d_list = tr3d._tr2d_list;
    if (tracer2d_list.size() != 3 || tracer2d_list[1]._pt_center != tr2d_4._pt_center)
    {
        std::cout << "test_function_2: updateTracer2D failed (line " << __LINE__ << ")" << std::endl;
        std::cout << "tracer2d_list.size() = " << tracer2d_list.size() << std::endl;
        return false;
    }

    tr3d.updateTracer2D({tr2d_4, tr2d_4}, {0, 2});
    tracer2d_list = tr3d._tr2d_list;
    if (tracer2d_list.size() != 2 || 
        tracer2d_list[0]._pt_center != tr2d_4._pt_center || 
        tracer2d_list[1]._pt_center != tr2d_4._pt_center)
    {
        std::cout << "test_function_2: updateTracer2D failed (line " << __LINE__ << ")" << std::endl;
        std::cout << "tracer2d_list.size() = " << tracer2d_list.size() << std::endl;
        return false;
    }
    ////////////////////////

    // test projectTracer2D
    std::vector<Camera> cam_list_all;
    std::vector<int> camid_list_2 = {0, 1, 2, 3};
    for (int i = 0; i < 5; i ++)
    {
        Camera cam ("../test/inputs/test_ObjectInfo/cam"+std::to_string(i+1)+".txt");
        cam_list_all.push_back(cam);
    }
    tr3d.projectTracer2D(camid_list_2, cam_list_all);
    tracer2d_list = tr3d._tr2d_list;
    camid_list = tr3d._camid_list;
    if (tracer2d_list.size() != 4 || camid_list.size() != 4)
    {
        std::cout << "test_function_2: projectTracer2D failed (line " << __LINE__ << ")" << std::endl;
        std::cout << "tracer2d_list.size() = " << tracer2d_list.size() << std::endl;
        std::cout << "camid_list.size() = " << camid_list.size() << std::endl;
        return false;
    }
    for (int i = 0; i < 4; i++)
    {
        if (tracer2d_list[i]._pt_center != cam_list_all[i].project(tr3d._pt_center) || camid_list[i] != camid_list_2[i])
        {
            std::cout << "test_function_2: projectTracer2D failed (line " << __LINE__ << ")" << std::endl;
            return false;
        }
    }
    ////////////////////////

    return true;
}


int main()
{
    IS_TRUE(test_function_1());
    IS_TRUE(test_function_2());

    return 0;
}