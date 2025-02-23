#include "ObjectInfo.h"

void Tracer3D::addTracer2D(Tracer2D const& tracer2d, int cam_id)
{
    _camid_list.push_back(cam_id);
    _tr2d_list.push_back(tracer2d);
    _n_2d ++;
}

void Tracer3D::addTracer2D(std::vector<Tracer2D> const& tracer2d_list, std::vector<int> const& camid_list)
{
    if (tracer2d_list.size() != camid_list.size())
    {
        std::cerr << "Tracer3D::addTracer2D: tracer2d_list.size() != camid_list.size()" << std::endl;
        throw error_size;
    }

    _camid_list.insert(_camid_list.end(), camid_list.begin(), camid_list.end());
    _tr2d_list.insert(_tr2d_list.end(), tracer2d_list.begin(), tracer2d_list.end());
    _n_2d += camid_list.size();
}

void Tracer3D::removeTracer2D(int cam_id)
{
    for (int i = 0; i < _n_2d; i ++)
    {
        if (_camid_list[i] == cam_id)
        {
            _camid_list.erase(_camid_list.begin() + i);
            _tr2d_list.erase(_tr2d_list.begin() + i);
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
    _tr2d_list.clear();
    _n_2d = 0;
}

void Tracer3D::updateTracer2D(Tracer2D const& tracer2d, int cam_id)
{
    for (int i = 0; i < _n_2d; i ++)
    {
        if (_camid_list[i] == cam_id)
        {
            _tr2d_list[i] = tracer2d;
            break;
        }
    }
}

void Tracer3D::updateTracer2D(std::vector<Tracer2D> const& tracer2d_list, std::vector<int> const& camid_list)
{
    if (tracer2d_list.size() != camid_list.size())
    {
        std::cerr << "Tracer3D::updateTracer2D: tracer2d_list.size() != camid_list.size()" << std::endl;
        throw error_size;
    }

    _tr2d_list = tracer2d_list;
    _camid_list = camid_list;
    _n_2d = camid_list.size();
}

void Tracer3D::projectObject2D(std::vector<int> const& camid_list, std::vector<Camera> const& cam_list_all)
{
    _n_2d = camid_list.size();
    _camid_list = camid_list;
    _tr2d_list.resize(_n_2d);

    int cam_id;
    for (int i = 0; i < _n_2d; i ++)
    {
        cam_id = _camid_list[i];
        _tr2d_list[i]._pt_center = cam_list_all[cam_id].project(_pt_center);
        _tr2d_list[i]._r_px = _r2d_px;
    }
}

void Tracer3D::getTracer2D(Tracer2D& tracer2d, int cam_id)
{
    for (int i = 0; i < _camid_list.size(); i ++)
    {
        if (_camid_list[i] == cam_id)
        {
            tracer2d = _tr2d_list[i];
            break;
        }
    }
}

void Tracer3D::saveObject3D(std::ofstream& output, int n_cam_all) const
{
    output << _pt_center[0] << "," << _pt_center[1] << "," << _pt_center[2] << "," << _error << "," << _n_2d;
    
    std::vector<double> pt2d_list(n_cam_all*2, IMGPTINIT);
    for (int i = 0; i < _n_2d; i ++)
    {
        pt2d_list[_camid_list[i]*2] = _tr2d_list[i]._pt_center[0];
        pt2d_list[_camid_list[i]*2+1] = _tr2d_list[i]._pt_center[1];
    }

    for (int i = 0; i < n_cam_all; i ++)
    {
        output << "," << pt2d_list[i*2] << "," << pt2d_list[i*2+1];
    }

    output << "\n";
}

