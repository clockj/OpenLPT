#include "ObjectInfo.h"
#include "SphereInfo.h"

void Sphere::addCircle(Circle const& circle, int cam_id)
{
    _camid_list.push_back(cam_id);
    _circle_list.push_back(circle);
    _n_2d ++;
}

void Sphere::addCircle(std::vector<Circle> const& circle_list, std::vector<int> const& camid_list)
{
    if (circle_list.size() != camid_list.size())
    {
        std::cerr << "Sphere::addCircle: cirlce_list.size() != camid_list.size()" << std::endl;
        throw error_size;
    }

    _camid_list.insert(_camid_list.end(), camid_list.begin(), camid_list.end());
    _circle_list.insert(_circle_list.end(), circle_list.begin(), circle_list.end());
    _n_2d += camid_list.size();
}

void Sphere::removeCircle(int cam_id)
{
    for (int i = 0; i < _n_2d; i ++)
    {
        if (_camid_list[i] == cam_id)
        {
            _camid_list.erase(_camid_list.begin() + i);
            _circle_list.erase(_circle_list.begin() + i);
            _n_2d --;
            break;
        }
    }
}

void Sphere::removeCircle(std::vector<int> const& camid_list)
{
    for (int i = 0; i < camid_list.size(); i ++)
    {
        removeCircle(camid_list[i]);
    }
}

void Sphere::clearCircle()
{
    _camid_list.clear();
    _circle_list.clear();
    _n_2d = 0;
}

void Sphere::updateCircle(Circle const& circle, int cam_id)
{
    for (int i = 0; i < _n_2d; i ++)
    {
        if (_camid_list[i] == cam_id)
        {
            _circle_list[i] = circle;
            break;
        }
    }
}

void Sphere::updateCircle(std::vector<Circle> const& circle_list, std::vector<int> const& camid_list)
{
    if (circle_list.size() != camid_list.size())
    {
        std::cerr << "Sphere::updateCircle: circle_list.size() != camid_list.size()" << std::endl;
        throw error_size;
    }

    _circle_list = circle_list;
    _camid_list = camid_list;
    _n_2d = camid_list.size();
}

void Sphere::projectObject2D(std::vector<int> const& camid_list, std::vector<Camera> const& cam_list_all)
{
    _n_2d = camid_list.size();
    _camid_list = camid_list;
    _circle_list.resize(_n_2d);

    int cam_id;
    for (int i = 0; i < _n_2d; i ++)
    {
        cam_id = _camid_list[i];
        _circle_list[i]._pt_center = cam_list_all[cam_id].project(_pt_center);
    }
}

void Sphere::getCircle(Circle& circle, int cam_id)
{
    for (int i = 0; i < _camid_list.size(); i ++)
    {
        if (_camid_list[i] == cam_id)
        {
            circle = _circle_list[i];
            break;
        }
    }
}

void Sphere::saveObject3D(std::ofstream& output, int n_cam_all) const
{
    // output: x, y, z, r_mm, error, n_2d, x1, y1, r1_px, x2, y2, r2_px, ...
    output << _pt_center[0] << "," << _pt_center[1] << "," << _pt_center[2] << ",";
    output << _r_mm << "," << _error << "," << _n_2d;
    
    std::vector<double> pt2d_list(n_cam_all*3, IMGPTINIT);
    for (int i = 0; i < _n_2d; i ++)
    {
        pt2d_list[_camid_list[i]*3]   = _circle_list[i]._pt_center[0];
        pt2d_list[_camid_list[i]*3+1] = _circle_list[i]._pt_center[1];
        pt2d_list[_camid_list[i]*3+2] = _circle_list[i]._r_px;
    }

    for (int i = 0; i < n_cam_all; i ++)
    {
        output << "," << pt2d_list[i*3] 
               << "," << pt2d_list[i*3+1] 
               << "," << pt2d_list[i*3+2];
    }

    output << "\n";
}
