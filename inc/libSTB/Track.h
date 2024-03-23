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


template<class T3D>
class Track
{
public:
    std::vector<T3D> _obj3d_list;
    std::vector<int> _t_list; // frame ID list
    int _n_obj3d = 0;
    bool _active = true;


    // Functions //
    Track() {};
    Track(T3D const& obj3d, int t);
    Track(Track const& track);
    ~Track() {};

    // Add a new point to the Track
    // only non-fake points should be added to the track
    void addNext(T3D const& obj, int t);

    // Add another Track onto the end of this one
    void addNext(Track const& t);
    
    // Predict the next position of the track
    void predictNext(T3D& obj3d);

    // write the track to a file
    void saveTrack(std::ofstream& output, int track_id, float fps = 1, int n_cam_all = 0);
    
    // member operators
    Track<T3D>& operator=(Track<T3D> const& t);

private:
    void predLMSWiener (T3D& obj3d);

};

#include "Track.hpp"

#endif