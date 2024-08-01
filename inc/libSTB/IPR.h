#ifndef IPR_H
#define IPR_H

#include <time.h>
#include <vector>
#include <algorithm>
#include <random>

#include "Matrix.h"
#include "Camera.h"
#include "ObjectInfo.h"
#include "ObjectFinder.h"
#include "StereoMatch.h"
#include "Shake.h"
#include "OTF.h"

struct IPRParam
{
    bool tri_only = false;
    int n_thread = 6; // number of threads
    int n_loop_ipr = 4; // 
    int n_loop_ipr_reduced = 2; // number of shake times after reducing cameras

    // Object finder parameters
    int n_obj2d_max = 1e5; // maximum number of tracers in each camera

    // Stereo match parameters
    double tol_2d = 1.; // [px]
    double tol_3d = 2.4e-2;  // [mm], 0.60 * (x_max-x_min)/1000
    int check_id = 3; // check_id <= n_use !
    double check_radius = 3; // [px]

    // Shake parameters
    int n_loop_shake = 1; // number of shake times, using gradient descent
    double shake_width = 2.4e-2; // [mm]
    double ghost_threshold = 0.1; // Ghost threshold: remove residue > mean + ghost_threshold * std
};


class IPR 
{
private:
    CamList& _cam_list;
    int _n_cam_all;

    void resetCamIDList ();

public:
    std::vector<Image> _imgRes_list;
    IPRParam _param;

    IPR (CamList& cam_list, std::vector<Image> const& imgOrig_list, IPRParam const& param) 
        : _cam_list(cam_list), _n_cam_all(cam_list.cam_list.size()), _imgRes_list(imgOrig_list), _param(param) 
    {};

    ~IPR () {};

    // Run IPR
    void runIPR(std::vector<Tracer3D>& tr3d_list_all, std::vector<double> const& tr2d_properties, OTF const& otf, int n_reduced = 0);

    // Run IPR with reduced cameras
    void reducedCamLoop (std::vector<Tracer3D>& tr3d_list_all, std::vector<double> const& tr2d_properties, 
    OTF const& otf, std::deque<std::vector<int>> const& cam_id_all, int n_cam);
    void createCamID (std::deque<std::vector<int>>& cam_id_all, std::vector<int> cam_id, int id, int n_rest);

    // Save objct info
    template <class T3D>
    void saveObjInfo (std::string const& filename, std::vector<T3D> const& obj3d_list)
    {
        std::ofstream file(filename);

        file << "WorldX,WorldY,WorldZ,Error,Ncam";
        for (int i = 0; i < _n_cam_all; i ++)
        {
            file << "," << "cam" << i << "_" << "x(col)" 
                 << "," << "cam" << i << "_" << "y(row)";
        }
        file << "\n";

        for (int i = 0; i < obj3d_list.size(); i ++)
        {
            obj3d_list[i].saveObject3D(file, _n_cam_all);
        }
        file.close();
    };
};

#include "IPR.hpp"

#endif // !IPR_H
