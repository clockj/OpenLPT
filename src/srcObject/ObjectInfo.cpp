#include "ObjectInfo.h"

void Tracer3D::addTracer2D(Tracer2D const& tracer2d, int cam_id)
{
    _camid_list.push_back(cam_id);
    _tracer2d_list.push_back(tracer2d);
    _n_2d ++;
}

void Tracer3D::addTracer2D(std::vector<Tracer2D> const& tracer2d_list, std::vector<int> const& camid_list)
{
    if (tracer2d_list.size() != camid_list.size())
    {
        std::cout << "Tracer3D::addTracer2D: tracer2d_list.size() != camid_list.size()" << std::endl;
        throw error_size;
    }

    _camid_list.insert(_camid_list.end(), camid_list.begin(), camid_list.end());
    _tracer2d_list.insert(_tracer2d_list.end(), tracer2d_list.begin(), tracer2d_list.end());
    _n_2d += camid_list.size();
}

void Tracer3D::removeTracer2D(int cam_id)
{
    for (int i = 0; i < _n_2d; i ++)
    {
        if (_camid_list[i] == cam_id)
        {
            _camid_list.erase(_camid_list.begin() + i);
            _tracer2d_list.erase(_tracer2d_list.begin() + i);
            _n_2d --;
            break;
        }
    }
}

void Tracer3D::removeTracer2D(std::vector<int> const& camid_list)
{
    for (int i = 0; i < camid_list.size(); i ++)
    {
        removeTracer2D(camid_list[i]);
    }
}

void Tracer3D::clearTracer2D()
{
    _camid_list.clear();
    _tracer2d_list.clear();
    _n_2d = 0;
}

void Tracer3D::updateTracer2D(Tracer2D const& tracer2d, int cam_id)
{
    for (int i = 0; i < _n_2d; i ++)
    {
        if (_camid_list[i] == cam_id)
        {
            _tracer2d_list[i] = tracer2d;
            break;
        }
    }
}

void Tracer3D::updateTracer2D(std::vector<Tracer2D> const& tracer2d_list, std::vector<int> const& camid_list)
{
    if (tracer2d_list.size() != camid_list.size())
    {
        std::cout << "Tracer3D::updateTracer2D: tracer2d_list.size() != camid_list.size()" << std::endl;
        throw error_size;
    }

    _tracer2d_list = tracer2d_list;
    _camid_list = camid_list;
    _n_2d = camid_list.size();
}

void Tracer3D::projectTracer2D(std::vector<int> const& camid_list, std::vector<Camera>& cam_list_all)
{
    _n_2d = camid_list.size();
    _camid_list = camid_list;
    _tracer2d_list.reserve(_n_2d);
    _tracer2d_list.resize(_n_2d);

    for (int i = 0; i < _n_2d; i ++)
    {
        int cam_id = _camid_list[i];
        _tracer2d_list[i]._pt_center = cam_list_all[cam_id].project(_pt_center);
    }
}

void Tracer3D::getTracer2D(Tracer2D& tracer2d, int cam_id)
{
    for (int i = 0; i < _camid_list.size(); i ++)
    {
        if (_camid_list[i] == cam_id)
        {
            tracer2d = _tracer2d_list[i];
            break;
        }
    }
}

void Tracer3D::getTracer2D(std::vector<Tracer2D>& tracer2d_list)
{
    tracer2d_list = _tracer2d_list;
}

void Tracer3D::getCamID(std::vector<int>& camid_list)
{
    camid_list = _camid_list;
}

