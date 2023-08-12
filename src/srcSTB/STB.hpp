#ifndef STB_HPP
#define STB_HPP

#include "STB.h"


template<class T>
STB<T>::STB (std::string file) :_ngrid_xyz(3,0)
{
    // Read param
    std::ifstream stb_config(file, std::ios::in);
    std::stringstream parsed;
    std::string line;
    while (std::getline(stb_config, line))
    {
        size_t commentpos = line.find('#');
        if (commentpos > 0)
        {
            if (commentpos < std::string::npos)
            {
                line.erase(commentpos);
            }
            parsed << line << '\t';
        }
    }
    stb_config.close();


    //        //
    // Camera //
    //        //
    parsed >> _n_cam;
    parsed >> _img_folder;
    for (int i = 0; i < _n_cam; i ++)
    {
        Camera cam(
            _img_folder 
            + "cam" + std::to_string(i+1) 
            + "Params.txt"
        );
        _cam_list.push_back(cam);
    }
    _n_pix_w = _cam_list[0].GetNpixw();
    _n_pix_h = _cam_list[1].GetNpixh();
    for (int i = 0; i < _n_cam; i ++)
    {
        ImageIO img;
        img.LoadImgPath(
            _img_folder, 
            "cam" + std::to_string(i+1) + "ImageNames.txt"
        );
        _imgio_list.push_back(img);
    }


    //          //
    // Tracking //
    //          //
    parsed >> _first; _first -= 1;
    parsed >> _last;  _last -= 1;

    // View area
    parsed >> _xyz_limit._x_min;
    parsed >> _xyz_limit._x_max;
    parsed >> _xyz_limit._y_min;
    parsed >> _xyz_limit._y_max;
    parsed >> _xyz_limit._z_min;
    parsed >> _xyz_limit._z_max;
    _otf = OTF(4, 3, _xyz_limit);

    parsed >> _vox_to_mm;

    // Initial phase
    parsed >> _ipr_flag; // 1 for using ipr in initialphase (0 for use .csv files)
    parsed >> _r_search_init; _r_search_init *= _vox_to_mm; // mm

    // Convergence phase
    parsed >> _shake_shift; _shake_shift *= _vox_to_mm;
    parsed >> _space_avg; _space_avg *= _vox_to_mm;
    parsed >> _shift_betframe_max; _shift_betframe_max *= _vox_to_mm;
    parsed >> _shift_abs_max; _shift_abs_max *= _vox_to_mm;
    parsed >> _shift_rel_max;
    parsed >> _fpt;
    parsed >> _int_lower;


    //                  //
    // Predictive Field //
    //                  //
    int n_grid;
    parsed >> n_grid; _ngrid_xyz[0] = n_grid;
    parsed >> n_grid; _ngrid_xyz[1] = n_grid;
    parsed >> n_grid; _ngrid_xyz[2] = n_grid;
    parsed >> _r_search_pred; _r_search_pred *= _vox_to_mm;

    //     //
    // IPR //
    //     //
    // Load IPR option
    int option;
    parsed >> option;
    if (option == 1)
    {
        _TRIONLY = true;
    }
    else if (option == 2)
    {
        _IPRONLY = true;
    }

    // Load IPR iteration
    parsed >> _ipr_loop_outer;
    parsed >> _ipr_loop_inner;
    parsed >> _ipr_int_lower;
    parsed >> _ipr_tol_2d; _ipr_tol_2d *= _vox_to_mm;
    parsed >> _ipr_tol_3d; _ipr_tol_3d *= _vox_to_mm;

    parsed >> _ipr_is_reduced;
    parsed >> _ipr_loop_outer_reduced;
    parsed >> _ipr_tol_2d_reduced; _ipr_tol_2d_reduced *= _vox_to_mm;
    parsed >> _ipr_tol_3d_reduced; _ipr_tol_3d_reduced *= _vox_to_mm;
    
    //             //
    // Object Info //
    //             //
    if (typeid(T) == typeid(TracerInfo))
    {   
        SetTracerParam (parsed);
    }
    else 
    {
        std::cerr << "StereoMatch: " 
                  << "class " << typeid(T).name()
                  << " is not included in StereoMatch!"
                  << std::endl;
        throw;
    }
}


template<class T>
void STB<T>::SetTracerParam (std::stringstream& parsed)
{
    // Load particle radius
    int r_tracer;
    parsed >> r_tracer;
    _obj_info_sample.SetRadiusPixel(r_tracer);

    _imgint_max_list = std::vector<int>(_n_cam, 0);
    _imgint_min_list = std::vector<int>(_n_cam, 0);
    for (int i = 0; i < _n_cam; i ++)
    {
        parsed >> _imgint_max_list[i];
    }
    for (int i = 0; i < _n_cam; i ++)
    {
        parsed >> _imgint_min_list[i];
    }
}


template<class T>
void STB<T>::InitialPhase ()
{
    int endframe = _last + 1;
    if (!_TRIONLY && !_IPRONLY)
    {
        endframe = (_last >= _first + 4) ? _first+4 : _last+1;
    }

    Matrix<double> intensity(_n_pix_h, _n_pix_w);
    for (int i = 0; i < _n_cam; i ++)
    {
        _img_org_list.push_back(intensity);
        // _img_res_list.push_back(intensity);
    }

    for (int frame = _first; frame < endframe; frame ++)
    {
        std::cout << "IPR on frame " << frame << " of " << endframe << std::endl;

        if (_ipr_flag)
        {
            // Load img
            for (int i = 0; i < _n_cam; i ++)
            {
                _img_org_list[i] = _imgio_list[i].LoadImg(frame);
            }

            // IPR
            IPR<T> ipr(_img_org_list, _imgint_max_list, _imgint_min_list, _cam_list, _otf);
            ipr.SetTRIONLY(_TRIONLY);
            ipr.SetIPRTimes(_ipr_loop_outer);
            ipr.SetShakeTimes(_ipr_loop_inner);
            ipr.SetTol2D(_ipr_tol_2d);
            ipr.SetTol3D(_ipr_tol_3d);

            std::vector<T> obj_list;
            ipr.RunIPR(obj_list, _ipr_is_reduced, 1);
            _ipr_matched.push_back(obj_list);
        }
    }

    if (!_TRIONLY && !_IPRONLY)
    {
        for (int frame = _first; frame < endframe-1; frame ++)
        {   
            int currframe = frame;
            int nextframe = frame + 1;

            std::cout << "STB initial phase tracking b/w frames: "
                      << currframe << " & " << nextframe << "; ";
            
            PredField<T> pf(_xyz_limit, _ngrid_xyz, _ipr_matched[currframe-_first], _ipr_matched[nextframe-_first], _r_search_pred);

            // link particles b/w two frames
            clock_t t_start, t_end;
            t_start = clock();

            // Extend all the tracks that are active in current frame
            Matrix<double> vel_curr(3,1);
            for (typename std::deque<Track<T>>::iterator tr = _active_short_track.begin(); tr != _active_short_track.end();)
            {
                bool active;
                vel_curr = pf.PtInterp(tr->Last().GetCenterPos());
                MakeLink(nextframe, vel_curr, _r_search_pred, *tr, active);

                if (!active)
                {
                    tr = _active_short_track.erase(tr);
                }
                else
                {
                    ++ tr;
                }
            }

            // Start a track for all particles left untracked in current frame
            StartTrack(currframe, pf);

            t_end = clock();
            std::cout << "; link particles time: " 
                      << (double) (t_end - t_start)/CLOCKS_PER_SEC
                      << std::endl;
        }

        // Move all tracks longer than 3 from _active_short_track to _active_long_track
        for (typename std::deque<Track<T>>::iterator tr = _active_short_track.begin(); tr != _active_short_track.end();)
        {
            if (tr->Length() > 3)
            {
                _active_long_track.push_back(*tr);
                tr = _active_short_track.erase(tr);
            }
            else
            {
                ++ tr;
            }
        }

        // Add the untracked particles from frame 4 to _active_short_track
        int m = (endframe-1)-_first;
        for (int i = 0; i < _ipr_matched[m].size(); i ++)
        {
            if (!_ipr_matched[m][i].IsTrack())
            {
                Track<T> temp(_ipr_matched[m][i], endframe-1);
                _active_short_track.push_back(temp);
            }
        }

        std::cout << "Done with initial phase!" << std::endl;
        std::cout << "\tNo. of active short tracks: " << _active_short_track.size() 
                  << "\tNo. of active long tracks: " << _active_long_track.size()
                  << "\tNo. of exited tracks: " << _exit_track.size()
                  << std::endl << std::endl;
    }
    
}


template<class T>
void STB<T>::Run ()
{
    InitialPhase();
}


template<class T>
void STB<T>::MakeLink(int nextframe, const Matrix<double>& vel_curr, double radius, Track<T>& track, bool& active)
{
    Matrix<double> pt_last(3,1);
    Matrix<double> pt_estimate(3,1);

    pt_last = track.Last().GetCenterPos();
    pt_estimate = pt_last + vel_curr;

    int len = track.Length();
    int obj_id = _UNLINKED;
    int m = nextframe-_first;

    NearestNeighbor(_ipr_matched[m], radius, pt_estimate, obj_id);

    if (obj_id == _UNLINKED)
    {
        active = false;
    }
    else
    {
        _ipr_matched[m][obj_id].SetIsTrack(true);
        track.AddNext(_ipr_matched[m][obj_id], nextframe);
        active = true;
    }
}


template<class T>
void STB<T>::NearestNeighbor(std::vector<T>& obj_list, double radius, const Matrix<double>& pt_estimate, int& obj_id)
{
    obj_id = -1;
    double min2 = radius * radius;
    double dis2 = 0;

    Matrix<double> pt(3,1);
    for (int i = 0; i < obj_list.size(); i ++)
    {
        pt = obj_list[i].GetCenterPos();
        // dis2 = pt_estimate.DistSqr(pt);
        dis2 = pt.DistSqr(pt_estimate);

        if (dis2 < min2)
        {
            obj_id = i;
            min2 = dis2;
        }
    }
}


template<class T>
void STB<T>::StartTrack(int frame, PredField<T>& pf)
{
    int m = frame-_first;
    Matrix<double> vel_curr(3,1);
    for (int i = 0; i < _ipr_matched[m].size(); i ++)
    {
        if (!_ipr_matched[m][i].IsTrack())
        {
            Track<T> init_tr;
            init_tr.AddNext(_ipr_matched[m][i], frame);
            _ipr_matched[m][i].SetIsTrack(true);

            vel_curr = pf.PtInterp(_ipr_matched[m][i].GetCenterPos());

            bool active;
            MakeLink(frame+1, vel_curr, _r_search_pred, init_tr, active);

            if (active)
            {
                _active_short_track.push_back(init_tr);
            }
        }
    }
}


#endif