#ifndef TRACK_H
#define TRACK_H

// #include <deque>
#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <typeinfo>

#include "STBCommons.h"
#include "Matrix.h"
#include "ObjectInfo.h"
#include "KalmanFilter.h"

// Predictor param
struct TrackPredParam
{
    TrackPredID type = WIENER;
    std::vector<double> param;    

    // Wiener filter: no param
    // Kalman filter: 
    //  _param = {
    //      sigma_q, # process noise
    //      sigma_rx, sigma_ry, sigma_rz, # measurement noise
    //      sigma_p0_x, sigma_p0_y, sigma_p0_z, 
    //      sigma_p0_vx, sigma_p0_vy, sigma_p0_vz # initial state covariance
    //  }
};

template<class T3D>
class Track
{
public:
    std::vector<T3D> _obj3d_list;
    std::vector<int> _t_list; // frame ID list
    int _n_obj3d = 0;
    bool _active = true;
    TrackPredParam _track_pred_param;

    // Functions //
    Track() {};
    Track(TrackPredParam const& track_pred_param);
    Track(T3D const& obj3d, int t);
    Track(T3D const& obj3d, int t, TrackPredParam const& track_pred_param);
    Track(Track const& track);
    ~Track() {};

    // Add a new point to the Track
    // only non-fake points should be added to the track
    void addNext(T3D const& obj, int t);

    // Add another Track onto the end of this one
    void addNext(Track const& t);
    
    // Predict the next position of the track
    // but not update the obj2d list
    void predictNext(T3D& obj3d);

    // Update the track with the new observation
    void update();

    // write the track to a file
    void saveTrack(std::ofstream& output, int track_id, float fps = 1, int n_cam_all = 0);
    
    // member operators
    Track<T3D>& operator=(Track<T3D> const& t);

private:
    void predLMSWiener (T3D& obj3d);

    // Kalman filter
    KalmanFilter _kf;
    bool _is_kf_init = false;
    void predKalman (T3D& obj3d);
};

#include "Track.hpp"

#endif