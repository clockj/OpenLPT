#ifndef SPHEREINFO_H
#define SPHEREINFO_H 

#include "ObjectInfo.h"

class Circle : public Object2D
{
public:
    double _r_px = 6; // [px], for shaking

    Circle () {};
    Circle (Circle const& circle) : Object2D(circle), _r_px(circle._r_px) {};
    Circle (Pt2D const& pt_center) : Object2D(pt_center) {};
    ~Circle () {};
};

// TODO: add a function to calculate _r_mm
// TODO: add a function to project _r_mm to _r_px 
class Sphere : public Object3D
{
public:
    int _n_2d = 0;
    double _error = 0; // [mm]
    double _r_mm = 0;  // [mm]
    std::vector<int> _camid_list; // min id = 0
    std::vector<Circle> _circle_list;

    Sphere () {};
    Sphere (Sphere const& sphere) : Object3D(sphere), _n_2d(sphere._n_2d), _error(sphere._error), _r_mm(sphere._r_mm), _camid_list(sphere._camid_list), _circle_list(sphere._circle_list) {};
    Sphere (Pt3D const& pt_center) : Object3D(pt_center) {};

    // user has to make sure the cam_id is unique
    void addCircle (Circle const& circle, int cam_id);
    void addCircle (std::vector<Circle> const& circle_list, std::vector<int> const& camid_list);

    void removeCircle (int cam_id);
    void removeCircle (std::vector<int> const& camid_list);
    void clearCircle ();
    
    void updateCircle (Circle const& circle, int cam_id);

    // _tr2d_list and _camid_list will be updated to tracer2d_list and camid_list
    void updateCircle (std::vector<Circle> const& circle_list, std::vector<int> const& camid_list);

    // project 3D tracer to 2D
    // It will project the 3D tracer to 2D for each camera in camid_list
    //  and store the 2D tracers in _tr2d_list.
    // input:
    //  camid_list: camera id list (set _camid_list = camid_list)
    //  cam_list_all: all camera parameters, camid = 0, 1, 2, ...
    void projectObject2D (std::vector<int> const& camid_list, std::vector<Camera> const& cam_list_all);

    void getCircle (Circle& circle, int cam_id);

    // save the 3D tracer to a file
    void saveObject3D (std::ofstream& output, int n_cam_all) const;
};

#endif