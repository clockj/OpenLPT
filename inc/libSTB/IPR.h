#ifndef IPR_H
#define IPR_H

#include <time.h>
#include <vector>

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

    // Stereo match parameters
    double tol_2d = 1.; // [mm], 0.25 * (x_max-x_min)/1000
    double tol_3d = 2.4e-2;  // [mm], 0.60 * (x_max-x_min)/1000
    int check_id = 3; // check_id <= n_use !
    double check_radius = 3; // [px]

    // Shake parameters
    int n_loop_shake = 1; // number of shake times, using gradient descent
    double shake_width = 2*2.4e-2; // [mm]
    double ghost_threshold = 3; // Ghost threshold: remove residue > mean + ghost_threshold * std
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
    void runIPR(std::vector<Tracer3D>& tr3d_list_all, std::vector<double> const& tr2d_properties, OTF const& otf, bool is_reduced = false, int n_reduced = 1);

    // Run IPR with reduced cameras
    void reducedCamLoop (std::vector<Tracer3D>& tr3d_list_all, std::vector<double> const& tr2d_properties, 
    OTF const& otf, std::vector<int> const& cam_id, int n_cam);
    void createCamID (std::deque<std::vector<int>>& cam_id_all, std::vector<int> cam_id, int id, int n_rest);

};


#endif // !IPR_H
