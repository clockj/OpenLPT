#ifndef TRACK_HPP
#define TRACK_HPP

#include "Track.h"

template<class T>
void Track<T>::AddNext(const T& p, int t)
{
    _obj_list.push_back(p);
    _t_list.push_back(t);
    _n_obj ++;
    _active = true;
}

template<class T>
void Track<T>::AddNext(const Track& t)
{
    _obj_list.insert(_obj_list.end(), t._obj_list.begin(), t._obj_list.end());
    _t_list.insert(_t_list.end(), t._t_list.begin(), t._t_list.end());
    _n_obj = _obj_list.size();
}

template<class T>
void Track<T>::DeleteBack() 
{
    _obj_list.pop_back();
    _t_list.pop_back();
    _n_obj --;
}

template<class T>
void Track<T>::DeleteFront() 
{
    _obj_list.pop_front();
    _t_list.pop_front();
    _n_obj --;
}

template<class T>
void Track<T>::AddFront(const T& pt, int t)
{
    _obj_list.push_front(pt);
    _t_list.push_front(t);
    _n_obj ++;
}

template<class T>
void Track<T>::AddFront(const Track& t)
{
    _obj_list.insert(_obj_list.begin(), t._obj_list.begin(), t._obj_list.end());
    _t_list.insert(_t_list.begin(), t._t_list.begin(), t._t_list.end());
    _n_obj = _obj_list.size();
}

template<class T>
int Track<T>::Length() const
{
    if (_obj_list.empty()) 
    {
        return 0;
    }

    // calculate effective size of the track: if the track ends with one or
    // more estimated positions, don't count them
    int adjust = 0;
    for (int i = _n_obj-1; i > -1; i --, adjust ++) 
    {
        if (!_obj_list[i].IsFake()) 
        {
            break;
        }
    }

    return (_n_obj - adjust);
}

template<class T>
int Track<T>::GetTime(int index) const 
{
    return _t_list[index];
}

template<class T>
Matrix<double> Track<T>::GetPos(int n) const 
{
    return _obj_list[n].GetCenterPos();
}

template<class T>
T Track<T>::GetObj(int n) const 
{
    return _obj_list[n];
}

template<class T>
int Track<T>::NumFake() const
{
    int count = 0;
    for (int i = 0; i < Length(); i++) 
    {
        if (_obj_list[i].IsFake()) 
        {
            count++;
        }
    }

    return count;
}

template<class T>
void Track<T>::WriteGDF(std::ofstream& output, float index, float fps /* = 1 */, int n_cam /*=4*/) const
{
    int len = Length();
    Matrix<double> pt(3,1);
    for (int i = 0; i < len; ++i) 
    {
        // Format:
        // Track Index
        // X
        // Y
        // Z
        // Framenumber
        // Camera 1-4 X and Y
        // Info
        // Fake bit
        pt = _obj_list[i].GetCenterPos();
        output.write(reinterpret_cast<const char*>(&index), 4);
        float tmp = pt(0,0);
        output.write(reinterpret_cast<const char*>(&tmp), 4);
        tmp = pt(1,0);
        output.write(reinterpret_cast<const char*>(&tmp), 4);
        tmp = pt(2,0);
        output.write(reinterpret_cast<const char*>(&tmp), 4);
        tmp = static_cast<float>(_t_list[i]) / fps;
        output.write(reinterpret_cast<const char*>(&tmp), 4);

        for (int cam_id = 0; cam_id < n_cam; cam_id ++)
        {
            pt = _obj_list[i].GetMatchPosInfo(cam_id);
            tmp = pt(0,0);
            output.write(reinterpret_cast<const char*>(&tmp), 4);
            tmp = pt(1,0);
            output.write(reinterpret_cast<const char*>(&tmp), 4);
        }

        // tmp = _obj_list[i].Info();
        // output.write(reinterpret_cast<const char*>(&tmp), 4);
        if (_obj_list[i].IsFake()) 
        {
            tmp = 1;
        } 
        else 
        {
            tmp = 0;
        }
        output.write(reinterpret_cast<const char*>(&tmp), 4);
    }
}

template<class T>
Track<T>& Track<T>::operator=(const Track<T>& t)
{
    _obj_list = t._obj_list;
    _t_list = t._t_list;
    _occluded = t._occluded;
    _n_obj = t._n_obj;
    _active = t._active;
    
    return *this;
}

template<class T>
bool Track<T>::Exists(int t) 
{
    bool exists = false;
    if (t <= _t_list[_n_obj - 1] && t >= _t_list[0])
    {
        exists = true;
    }
    return exists;
}

template<class T>
void Track<T>::Clear()
{
    _obj_list.clear();
    _t_list.clear();
    _n_obj = 0;
    _occluded = 0;
}

#endif