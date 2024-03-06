#include "test.h"
#include "Camera.h"
#include "ImageIO.h"
#include "ObjectInfo.h"
#include "ObjectFinder.h"
#include "StereoMatch.h"

// test createIDMap
bool test_function_1 ()
{    
    CamList cam_list;
    for (int i = 0; i < 4; i ++)
    {
        Camera cam("../test/inputs/test_StereoMatch/cam" + std::to_string(i+1) + ".txt");
        cam_list.cam_list.push_back(cam);
        cam_list.useid_list.push_back(i);
    }

    // load image path
    std::vector<ImageIO> imgio_list;
    for (int i = 0; i < 4; i ++)
    {
        ImageIO imgio;
        imgio.loadImgPath("../test/inputs/test_StereoMatch/", "cam" + std::to_string(i+1) + "ImageNames" + ".txt");
        imgio_list.push_back(imgio);
    }

    std::cout << "here 1" << std::endl;

    // load image
    std::vector<Image> img_list;
    for (int i = 0; i < 4; i ++)
    {
        img_list.push_back(imgio_list[i].loadImg(0));
    }

    std::cout << "here 2" << std::endl;

    // find 2d tracer
    std::vector<std::vector<Tracer2D>> tr2d_list_all;
    
    std::vector<double> properties = {255, 30};
    ObjectFinder2D objfinder;

    for (int i = 0; i < 4; i ++)
    {
        std::vector<Tracer2D> tr2d_list;
        objfinder.findObject2D(tr2d_list, img_list[i], properties);
        tr2d_list_all.push_back(tr2d_list);
    }

    std::cout << "here 3" << std::endl;

    StereoMatchParam param;
    param.search_radius = 10;
    param.tor_2d = 0.5;
    param.tor_3d = 2.4e-2;
    param.n_thread = 4;
    param.fov_size = 40;
    param.is_delete_ghost = false;
    StereoMatch stereo_match(param, cam_list);

    std::cout << "here 4" << std::endl;

    stereo_match.createObjIDMap(tr2d_list_all);

    std::cout << "here 5" << std::endl;

    return true;
}


int main ()
{
    IS_TRUE(test_function_1());

    return 0;
}