#ifndef TRACK_HPP
#define TRACK_HPP

#include "Track.h"

template<class T3D>
Track<T3D>::Track(TrackPredParam const& track_pred_param) : _track_pred_param(track_pred_param) {} 

template<class T3D>
Track<T3D>::Track(T3D const& obj3d, int t) : _n_obj3d(1) 
{ 
    _obj3d_list.push_back(obj3d); 
    _t_list.push_back(t); 
}

template<class T3D>
Track<T3D>::Track(T3D const& obj3d, int t, TrackPredParam const& track_pred_param) : _n_obj3d(1), _track_pred_param(track_pred_param)
{ 
    _obj3d_list.push_back(obj3d); 
    _t_list.push_back(t); 
}

template<class T3D>
Track<T3D>::Track(Track const& track) 
    : _obj3d_list(track._obj3d_list), _t_list(track._t_list), _n_obj3d(track._n_obj3d), _active(track._active), _track_pred_param(track._track_pred_param), _kf(track._kf), _is_kf_init(track._is_kf_init) {}


template<class T3D>
void Track<T3D>::addNext(T3D const& obj, int t)
{
    _obj3d_list.push_back(obj);
    _t_list.push_back(t);
    _n_obj3d ++;
}

template<class T3D>
void Track<T3D>::addNext(Track const& track)
{
    _obj3d_list.insert(_obj3d_list.end(), track._obj3d_list.begin(), track._obj3d_list.end());
    _t_list.insert(_t_list.end(), track._t_list.begin(), track._t_list.end());
    _n_obj3d = _obj3d_list.size();
}

template<class T3D>
void Track<T3D>::saveTrack(std::ofstream& output, int track_id, float fps, int n_cam_all)
{    
    if (_n_obj3d != _obj3d_list.size())
    {
        std::cerr << "Track<Tracer3D>::saveTrack error at line " << __LINE__ << ": _n_obj3d != _obj3d_list.size()" << std::endl;
        std::cerr << "track_id: " << track_id << std::endl;
        std::cerr << "_n_obj3d: " << _n_obj3d << std::endl;
        std::cerr << "_obj3d_list.size(): " << _obj3d_list.size() << std::endl;
        return;
    }

    for (int i = 0; i < _n_obj3d; i ++)
    {
        output << track_id << "," << _t_list[i]/fps << ",";

        // print 2d info
        _obj3d_list[i].saveObject3D(output, n_cam_all);
    }
}

template<class T3D>
Track<T3D>& Track<T3D>::operator=(Track<T3D> const& track)
{
    _obj3d_list = track._obj3d_list;
    _t_list = track._t_list;
    _n_obj3d = track._n_obj3d;
    _active = track._active;
    _track_pred_param = track._track_pred_param;
    _kf = track._kf;
    _is_kf_init = track._is_kf_init;
    return *this;
}

template<class T3D>
void Track<T3D>::predictNext(T3D& obj3d)
{
    if (typeid(T3D) == typeid(Tracer3D))
    {
        if (_track_pred_param.type == WIENER)
        {
            predLMSWiener(obj3d);
        }
        else if (_track_pred_param.type == KALMAN)
        {
            predKalman(obj3d);
        }
        else
        {
            std::cerr << "Error: predictNext not implemented for this type" << std::endl;
        }
    }
    else
    {
        std::cerr << "Error: predictNext not implemented for this type" << std::endl;
    }
}

template<class T3D>
void Track<T3D>::predLMSWiener(T3D& obj3d)
{
    int order = 0;

    if (_n_obj3d < 3)
    {
        std::cerr << "Track<T3D>::predLMSWiener error at line " << __LINE__ << ": no enough points to predict " << _n_obj3d << "<" << 3 << std::endl;
        return;
    }
    // else if (_n_obj3d < 4)
    // {
    //     order = 3;
    // }
    else if (_n_obj3d < 6)
    {
        order = _n_obj3d - 1;
    }
    else 
    {
        order = 5;
    }

    // predict at each direction
    std::vector<double> series (order+1, 0);
    std::vector<double> filter_param (order, 0);
    for (int i = 0; i < 3; i ++)
    {
        for (int j = 0; j < order+1; j ++)
        {
            series[j] = _obj3d_list[_n_obj3d-1-order + j]._pt_center[i];
        }
        
        // Wiener filter does badly near zero
        // Make a shift to avoid zero-plane prediction
        bool shift_label = false;
        double shift = 10;
        if (std::fabs(series[order]) < 1)
        {
            shift_label = true;
            for (int j = 0; j < order+1; j ++)
            {
                series[j] += shift;
            }
        }
        
        // calculate the step
        double sum = 0;
        for (int j = 0; j < order; j ++)
        {
            sum += (series[j] * series[j]);
        }
        double step = 1 / sum;

        // initialize the filter param
        std::fill(filter_param.begin(), filter_param.end(), 0);
        double prediction = 0;
        double error = series[order] - prediction;

        int iter = 0;
        while (std::fabs(error) > SMALLNUMBER && iter < WIENER_MAX_ITER)
        // while (std::fabs(error) > SMALLNUMBER)
        {
            // Update parameters and calculate the prediction using the new filter param
            prediction = 0;
            for (int j = 0; j < order; j ++)
            {
                filter_param[j] += (step * series[j] * error);
                prediction += (filter_param[j] * series[j]);
            }

            error = series[order] - prediction;
            iter ++;
        }
        // std::cout << "Prediction iter: " << iter << std::endl;

        prediction = 0;
        for (int j = 0; j < order; j ++)
        {
            prediction += (filter_param[j] * series[j+1]);
        }
        if (shift_label)
        {
            prediction -= shift;
        }
        obj3d._pt_center[i] = prediction;
    } 

    if (typeid(T3D) == typeid(Tracer3D))
    {
        obj3d._r2d_px = _obj3d_list[_n_obj3d-1]._r2d_px;
    } 
    else
    {
        std::cerr << "Error: predictNext not implemented for this type" << std::endl;
    }

}

template<class T3D>
void Track<T3D>::update()
{
    if (_n_obj3d < 2)
    {
        std::cerr << "Track<T3D>::update error at line " << __LINE__ << ": no enough points to update " << _n_obj3d << "<" << 2 << std::endl;
        return;
    }

    if (typeid(T3D) == typeid(Tracer3D))
    {
        _kf.update(_obj3d_list[_n_obj3d-1]._pt_center);

        // // update the tracer position
        // for (int i = 0; i < 3; i ++)
        // {
        //     _obj3d_list[_n_obj3d-1]._pt_center[i] = _kf._x[i];
        // }
    }
    else
    {
        std::cerr << "Error: update not implemented for this type" << std::endl;
    }
}

template<class T3D>
void Track<T3D>::predKalman(T3D& obj3d)
{
    if (_n_obj3d < 2)
    {
        std::cerr << "Track<T3D>::predKalman error at line " << __LINE__ << ": no enough points to predict " << _n_obj3d << "<" << 2 << std::endl;
        return;
    }

    // initialize the Kalman filter
    if (!_is_kf_init)
    {
        _kf._F = Matrix<double> (6,6,0);
        for (int i = 0; i < 3; i ++)
        {
            _kf._F(i,i) = 1;
            _kf._F(i,i+3) = 1;
        }
        for (int i = 3; i < 6; i ++)
        {
            _kf._F(i,i) = 1;
        }
        
        _kf._H = Matrix<double> (3,6,0);
        for (int i = 0; i < 3; i ++)
        {
            _kf._H(i,i) = 1;
        }

        int param_id = 0;
        double sigma = _track_pred_param.param[param_id];
        _kf._Q = Matrix<double> (6,6,0);
        for (int i = 0; i < 3; i ++)
        {
            _kf._Q(i,i) = 0.25;
            _kf._Q(i,i+3) = 0.5;
        }
        for (int i = 3; i < 6; i ++)
        {
            _kf._Q(i,i) = 1;
            _kf._Q(i,i-3) = 0.5;
        }
        _kf._Q *= (sigma*sigma);
        param_id ++;

        _kf._R = Matrix<double> (3,3,0);
        for (int i = 0; i < 3; i ++)
        {
            sigma = _track_pred_param.param[param_id];
            _kf._R(i,i) = sigma * sigma;
            param_id ++;
        }

        _kf._P = Matrix<double> (6,6,0);
        for (int i = 0; i < 6; i ++)
        {
            sigma = _track_pred_param.param[param_id];
            _kf._P(i,i) = sigma * sigma;
            param_id ++;
        }

        _kf._x = Matrix<double> (6,1,0);
        for (int i = 0; i < 3; i ++)
        {
            _kf._x(i,0) = _obj3d_list[1]._pt_center[i];
        }
        for (int i = 3; i < 6; i ++)
        {
            _kf._x(i,0) = _obj3d_list[1]._pt_center[i-3] - _obj3d_list[0]._pt_center[i-3];
        }

        // warm up the filter
        for (int i = 1; i < _n_obj3d-1; i ++)
        {
            _kf.predict();
            _kf.update(_obj3d_list[i+1]._pt_center);
        }

        _is_kf_init = true;
    }
    else
    {
        update();
    }

    _kf.predict();

    // predict the tracer position
    for (int i = 0; i < 3; i ++)
    {
        obj3d._pt_center[i] = _kf._x[i];
    }
}


#endif