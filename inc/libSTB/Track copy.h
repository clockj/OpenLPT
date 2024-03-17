#ifndef TRACK_H
#define TRACK_H

// #include <deque>
#include <vector>
#include <iostream>
#include <fstream>

#include "STBCommons.h"
#include "Matrix.h"
#include "ObjectInfo.h"


template<class T3D>
class Track
{
private:
    std::vector<T3D> _obj3d_list;
    std::vector<int> _t_list;	// the time (as an integer frame number)
    int _occluded = 0; // a counter keeping track of the number of frames this 
                   // track hasn't had a real particle added to it.
    int _n_obj = 0;
    bool _active = true;

public:
    Track() {};
    Track(const T3D& obj3d, int t) : _n_obj(1) { _obj_list.push_back(obj3d); _t_list.push_back(t); };
    Track(const Track& track) : _obj_list(t._obj_list), _t_list(t._t_list), _occluded(t._occluded), _n_obj(t._n_obj), _active(t._active) {};

    ~Track() {};

    // Add a new point to the Track
    void AddNext(const T3D& obj, int t);
    // Add another Track onto the end of this one
    void AddNext(const Track& t);

    // Delete a point to the end of the Track
    void DeleteBack();

    void DeleteFront();

    // Add a new point in the begining Track
    void AddFront(const T3D& obj, int t);
    // Add another Track onto the begining of this one
    void AddFront(const Track& t);

    const T3D First() const {return _obj_list[0];}; // be careful: no checking for range here
    const T3D Second() const {return _obj_list[1];}; // be careful: no checking for range here
    const T3D Third() const {return _obj_list[2];}; // be careful: no checking for range here

    // return the last point on the Track
    const T3D Last() const {return _obj_list[_n_obj-1];}; // be careful: no checking for range here
    // return the penultimate point on the Track
    const T3D Penultimate() const {return _obj_list[_n_obj-2];}; // be careful: no checking for range here
    // return the antepenultimate point on the Track
    const T3D Antepenultimate() const {return _obj_list[_n_obj-3];}; // be careful: no checking for range here
    

    // get the length of the Track (not including any ending extrapolated points)
    int Length() const;
    // get the time of a particular element
    int GetTime(int index) const;
    // get the position of nth element
    Matrix<double> GetPos(int n) const;
    T3D GetObj(int n) const;

    // check the occlusion counter
    int OcclusionCount() const {return _occluded;};
    // increment the occlusion counter
    void Occluded() {_occluded ++;};
    // reset the occlusion counter
    void ResetCounter() {_occluded = 0;};

    // find the number of usable fake points on the track
    int NumFake() const;

    // write the track as part of a GDF file (see Tracker.h for format)
    void WriteGDF(std::ofstream& output, float index, float fps = 1, int n_cam = 4) const;
    
    // member operators
    Track<T3D>& operator=(const Track<T3D>& t);
    // return the nth point on the Track
    const T3D operator[](int n) {return _obj_list[n];};

    // print only the estimated points
    void PrintEstimates(std::ostream& os) const;
    
    // set track as active / inactive
    void SetActive(bool activity) { _active = activity; };
    
    // check if the track is active / inactive
    bool IsActive() { return _active; };

    // check if the track has a particle at time 't'
    bool Exists(int t);

    // free memory
    void Clear();

};

#include "Track.hpp"

#endif