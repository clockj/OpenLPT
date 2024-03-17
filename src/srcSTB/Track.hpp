#ifndef TRACK_HPP
#define TRACK_HPP

#include "Track.h"

template<class T3D>
Track<T3D>::Track(const T3D& obj3d, int t) : _n_obj3d(1) 
{ 
    _obj3d_list.push_back(obj3d); 
    _t_list.push_back(t); 
}

template<class T3D>
Track<T3D>::Track(const Track& track) 
    : _obj3d_list(track._obj3d_list), _t_list(track._t_list), _occluded(track._occluded), _n_obj3d(track._n_obj3d), _active(track._active) {}


template<class T3D>
void Track<T3D>::addNext(const T3D& obj, int t)
{
    _obj3d_list.push_back(obj);
    _t_list.push_back(t);
    _n_obj3d ++;
}

template<class T3D>
void Track<T3D>::addNext(const Track& track)
{
    _obj3d_list.insert(_obj3d_list.end(), track._obj3d_list.begin(), track._obj3d_list.end());
    _t_list.insert(_t_list.end(), track._t_list.begin(), track._t_list.end());
    _n_obj3d = _obj3d_list.size();
}

template<>
void Track<Tracer3D>::saveTrack(std::ofstream& output, int track_id, float fps, int n_cam_all)
{
    std::vector<double> pt2d_list(n_cam_all*2, -10);
    
    int cam_id;
    for (int i = 0; i < _n_obj3d; i ++)
    {
        output << track_id << "," << _t_list[i]/fps << "," << _obj3d_list[i]._pt_center[0] << "," << _obj3d_list[i]._pt_center[1] << "," << _obj3d_list[i]._pt_center[2] << "," << _obj3d_list[i]._error;

        // print 2d info
        std::fill(pt2d_list.begin(), pt2d_list.end(), -10);
        for (int j = 0; j < _obj3d_list[i]._n_2d; j ++)
        {
            cam_id = _obj3d_list[i]._camid_list[j];

            pt2d_list[cam_id*2] = _obj3d_list[i]._tr2d_list[j]._pt_center[0];
            pt2d_list[cam_id*2+1] = _obj3d_list[i]._tr2d_list[j]._pt_center[1];
        }

        for (int j = 0; j < n_cam_all; j ++)
        {
            output << "," << pt2d_list[j*2] << "," << pt2d_list[j*2+1];
        }
        output << "\n";
    }
}

template<class T3D>
Track<T3D>& Track<T3D>::operator=(const Track<T3D>& track)
{
    _obj3d_list = track._obj3d_list;
    _t_list = track._t_list;
    _occluded = track._occluded;
    _n_obj3d = track._n_obj3d;
    _active = track._active;
    return *this;
}


#endif