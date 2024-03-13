//
//	ObjectInfo.h 
//
//	Base class for various image formats
//	All image format types should inherit from Image 
//
//	Created by Shijie Zhong 07/11/2022 
//

#ifndef OBJECTINFO_H
#define OBJECTINFO_H 

#include <vector>
#include "Matrix.h"
#include "Camera.h"
#include "STBCommons.h"

// 2D object classes
class Object2D
{
public:
    Pt2D _pt_center; 

    Object2D () {};
    Object2D (Object2D const& object) : _pt_center(object._pt_center) {};
    Object2D (Pt2D const& pt_center) : _pt_center(pt_center) {};
    virtual ~Object2D () {};
};

class Tracer2D : public Object2D
{
public:
    int _r_px = 2; // [px], for shaking

    Tracer2D () {};
    Tracer2D (Tracer2D const& tracer) : Object2D(tracer), _r_px(tracer._r_px) {};
    Tracer2D (Pt2D const& pt_center) : Object2D(pt_center) {};
    ~Tracer2D () {};
};


// 3D object classes
class Object3D
{
public:
    Pt3D _pt_center; 

    Object3D () {};
    Object3D (Object3D const& object) : _pt_center(object._pt_center) {};
    Object3D (Pt3D const& pt_center) : _pt_center(pt_center) {};
    virtual ~Object3D () {};
};

class Tracer3D : public Object3D
{
public:
    int _n_2d = 0;
    double _error = 0; // [mm]
    std::vector<int> _camid_list; // min id = 0
    std::vector<Tracer2D> _tr2d_list;

    Tracer3D () {};
    Tracer3D (Tracer3D const& tracer3d) : Object3D(tracer3d), _n_2d(tracer3d._n_2d), _error(tracer3d._error), _camid_list(tracer3d._camid_list), _tr2d_list(tracer3d._tr2d_list) {};
    Tracer3D (Pt3D const& pt_center) : Object3D(pt_center) {};
    ~Tracer3D () {};

    // user has to make sure the cam_id is unique
    void addTracer2D (Tracer2D const& tracer2d, int cam_id);
    void addTracer2D (std::vector<Tracer2D> const& tracer2d_list, std::vector<int> const& camid_list);

    void removeTracer2D (int cam_id);
    void removeTracer2D (std::vector<int> const& camid_list);
    void clearTracer2D ();
    
    void updateTracer2D (Tracer2D const& tracer2d, int cam_id);

    // _tr2d_list and _camid_list will be updated to tracer2d_list and camid_list
    void updateTracer2D (std::vector<Tracer2D> const& tracer2d_list, std::vector<int> const& camid_list);

    // project 3D tracer to 2D
    // It will project the 3D tracer to 2D for each camera in camid_list
    //  and store the 2D tracers in _tr2d_list.
    // input:
    //  camid_list: camera id list (set _camid_list = camid_list)
    //  cam_list_all: all camera parameters, camid = 0, 1, 2, ...
    void projectTracer2D (std::vector<int> const& camid_list, std::vector<Camera> const& cam_list_all);

    void getTracer2D (Tracer2D& tracer2d, int cam_id);
};

#endif
