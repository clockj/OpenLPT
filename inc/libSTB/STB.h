#ifndef STB_H
#define STB_H

#include "Matrix.h"
#include "Camera.h"
#include "ObjectInfo.h"
#include "ImageIO.h"
#include "ObjectFinder.h"
#include "StereoMatch.h"
#include "OTF.h"
#include "Shake.h"
#include "IPR.h"
#include "PredField.h"
#include "Track.h"

#include <vector>
#include <deque>
#include <typeinfo>
#include <omp.h>
#include <time.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <filesystem>

namespace fs = std::filesystem;


template<class T3D>
class STB
{
public:
    std::vector<std::vector<T3D>> _ipr_matched;

    std::deque<Track<T3D>> _short_track_active;
    std::deque<Track<T3D>> _long_track_active;
    std::deque<Track<T3D>> _long_track_inactive; // only save very long tracks
    std::deque<Track<T3D>> _exit_track;

    // FUNCTIONS //

    // constructor
    STB(int frame_start, int frame_end, float fps, double vx_to_mm, int n_thread, std::string const& output_folder, CamList const& cam_list, AxisLimit const& axis_limit,  std::string const& file);

    // destructor
    ~STB() {};

    enum TrackType{ Inactive = 0, ActiveShort = 1, ActiveLong = 2, Exit = 3, InactiveLong = 4, Buffer = 5};

    // Process STB on frame frame_id
    // return img_list: residue images
    void processFrame(int frame_id, std::vector<Image>& img_list, bool is_update_img = false);

    

    // save tracks of one status
    void saveTracks (std::string const& file, std::deque<Track<T3D>>& tracks);

    // save all tracks at any status
    void saveTracksAll(std::string const& folder, int frame);

private:
    int _first = 0; // first frame ID
    int _last = 0; // last frame ID
    float _fps = 1; // frame rate
    double _vx_to_mm; // voxel to mm, 1[vx]=_vx_to_mm[mm]
    int _n_thread;
    std::string _output_folder;
    
    CamList _cam_list;
    int _n_cam_all=0;
    AxisLimit _axis_limit;

    // Tracking parameters
    bool _ipr_flag;
    double _r_objSearch; // mm, find object around a short track between two frames
    double _r_trackSearch; // mm, find tracks around a track
    int _n_initPhase; // number of frames for initial phase

    // Shake parameters
    double _shake_width; // mm

    // Predict Field parameters
    PFParam _pf_param;

    // IPR parameters
    bool _ipr_only = false;
    IPRParam _ipr_param;
    int _n_reduced = 0;

    // Object parameters
    std::vector<double> _obj_param;
    OTF _otf;

    // Dummy variables to identify the no. of tracks added and subtracted 
    int _a_sa = 0, _a_la = 0, _s_sa = 0, _s_la = 0, _a_li = 0;


    // FUNCTIONS //
    void createFolder (std::string const& folder);

    void loadObjParam (std::stringstream& config);

    void runInitPhase (int frame_id, std::vector<Image>& img_list, bool is_update_img = false);

    void runConvPhase (int frame_id, std::vector<Image>& img_list, bool is_update_img = false);

    // make all possible links for tracks
    int makeLink (Track<T3D> const& track, int nextframe, Pt3D const& vel_curr, double radius);

    void startTrack (int frame, PredField& pf);

    int linkShortTrack (Track<T3D> const& track, std::vector<T3D> const& obj3d_list, int n_iter);

    bool checkLinearFit (Track<T3D> const& track);

    // find nearest neighbor around a position
    int findNN (std::vector<T3D> const& obj3d_list, Pt3D const& pt3d_est, double radius);

};

#include "STB.hpp"

#endif