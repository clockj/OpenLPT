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

    Matrix<double> intensity(_n_pix_h, _n_pix_w); // TODO: 1024*1024
    for (int i = 0; i < _n_cam; i ++)
    {
        _img_org_list.push_back(intensity);
        // _img_res_list.push_back(intensity);
    }

    for (int frame = _first; frame < endframe; frame ++)
    {
        std::cout << "\nIPR on frame " << frame << " of " << endframe << std::endl;

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
            ipr.SetShakeWidth(_shake_shift);

            std::vector<T> obj_list;
            _ipr_matched.push_back(obj_list);
            ipr.RunIPR(_ipr_matched.back(), _ipr_is_reduced, 1);
        }
    }

    if (!_TRIONLY && !_IPRONLY)
    {
        std::cout << std::endl;
        for (int frame = _first; frame < endframe-1; frame ++)
        {   
            int currframe = frame;
            int nextframe = frame + 1;

            std::cout << "STB initial phase tracking b/w frames: "
                      << currframe << " & " << nextframe << "; ";
            
            PredField<T> pf(_xyz_limit, _ngrid_xyz, _ipr_matched[currframe-_first], _ipr_matched[nextframe-_first], _r_search_pred);
            pf.SaveField(_img_folder + "PredField_" + std::to_string(currframe) + "&" + std::to_string(nextframe) + ".csv");

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
            std::cout << "; total link particles time: " 
                      << (double) (t_end - t_start)/CLOCKS_PER_SEC
                      << std::endl;
        }

        // Move all tracks longer than 3 from _active_short_track to _active_long_track
        for (typename std::deque<Track<T>>::iterator tr = _active_short_track.begin(); tr != _active_short_track.end();)
        {
            if (tr->Length() >= 4)
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
                  << "; No. of active long tracks: " << _active_long_track.size()
                  << "; No. of exited tracks: " << _exit_track.size()
                  << std::endl;
    }
    _ipr_matched.clear();


    // Save all the data
    std::string address = _img_folder + "Track/InitialTrack/";
    std::string s = std::to_string(endframe-1);
    std::string file;

    // save the active long tracks
    file = "ActiveLongTrack" + s + ".csv";
    SaveTrack(_active_long_track, address + file);

    // save the active short tracks
    file = "ActiveShortTrack" + s + ".csv";
    SaveTrack(_active_short_track, address + file);

    // //  save the inactive tracks
    // file = "InactiveTrack" + s + ".csv";
    // SaveTrack(_inactive_track, address + file);

    //  save the exit tracks
    file = "ExitTrack" + s + ".csv";
    SaveTrack(_exit_track, address + file);

    // save the inactive long tracks
    file = "InactiveLongTrack" + s + ".csv";
    SaveTrack(_inactive_long_track, address + file);
}


template<class T>
void STB<T>::ConvergencePhase ()
{
    int endframe = _last;
    std::string address = _img_folder + "Tracks/ConvergedTracks/";
    // Track frame by frame using Wiener / polynomial predictor along with shaking
    for (int currframe = _first+3; currframe < endframe-1; currframe ++)
    {
        // Initialize some variables
        int c1 = _active_short_track.size();
        int c2 = _active_long_track.size();
        int c3 = _inactive_track.size();
        int c4 = _inactive_long_track.size();
        _a_as = 0; _a_al = 0; _a_is = 0; _s_as1 = 0; _s_as2 = 0; _s_as3 = 0; _s_as4 = 0; _s_al = 0; _a_il = 0;

        int nextframe = currframe + 1;
        std::cout << "\nSTB convergence phase tracking at frame: " << nextframe << std::endl;


        // Prediction for active long tracks
        std::cout << "\tPrediction: ";
        clock_t t_start, t_end;
        t_start = clock();

        std::vector<Matrix<double>> est_pos;
        Prediction(nextframe, est_pos); 
        std::vector<int> est_int(est_pos.size(), 1);

        t_end = clock();
        std::cout << (t_end-t_start)/CLOCKS_PER_SEC << "s" << std::endl;
        std::cout << "\tDone!" << std::endl;


        // Load img
        for (int i = 0; i < _n_cam; i ++)
        {
            _img_org_list[i] = _imgio_list[i].LoadImg(nextframe);
        }


        // Shaking prediction
        int n_obj = est_pos.size();
        if (n_obj)
        {
            // Create obj_list based on est_pos
            std::vector<T> obj_info_match_list(n_obj);
            T obj;
            for (int i = 0; i < n_obj; i ++)
            {
                obj.SetCenterPos(est_pos[i]);
                obj.SetMatchPosInfo(_cam_list);
                obj_info_match_list[i] = obj;
            }

            // Shaking prediction
            // correcting the predicted positions by shaking, removing wrong / ambiguous predictions and updating the residual images
            t_start = clock();

            Shake<T> s(
                _img_org_list,
                obj_info_match_list, // also used as output
                _shake_shift, // unit: mm
                _cam_list,
                _otf
            );
            s.SetShakeTimes (_ipr_loop_inner);
            s.RunShake(_TRIONLY);
            s.GetResImg(_img_org_list);

            t_end = clock();
            std::cout << "\tShaking time: " << (t_end-t_start)/CLOCKS_PER_SEC << "s" << std::endl;

            // Add the corrected particle position in nextFrame to its respective track
            std::vector<int> tr_id_keep_back = s.GetKeepID_BackWards();
            
            int n_obj_shake = tr_id_keep_back.size();
            if (n_obj_shake != obj_info_match_list.size())
            {
                std::cout << "STB::ConvergencePhase issues:\n"
                          << "n_obj_shake=" << n_obj_shake << "\n"
                          << "obj_info_match_list.size=" << obj_info_match_list.size() << std::endl;
                throw error_size;
            }          
            for (int i = 0; i < n_obj_shake; i ++)
            {
                int id = tr_id_keep_back[n_obj_shake-i-1];
                _active_long_track[id].AddNext(obj_info_match_list[i], nextframe);
            }

            // Update track list status
            //  from back to front
            std::vector<int> tr_id_rem_back = s.GetRemID();
            for (int i = 0; i < tr_id_rem_back.size(); i ++)
            {
                int id = tr_id_rem_back[i];
                if (_active_long_track[id].Length() > 6)
                {
                    _inactive_long_track.push_back(_active_long_track[id]);
                    _a_il ++;
                }
                else
                {               
                    _a_is ++;
                }
                _active_long_track.erase(_active_long_track.begin()+id);
                _s_al ++;
            }

            std::cout << "\tDone!" << std::endl;
        }


        // IPR on residues
        IPR<T> ipr(_img_org_list, _imgint_max_list, _imgint_min_list, _cam_list, _otf);
        ipr.SetTRIONLY(_TRIONLY);
        ipr.SetIPRTimes(_ipr_loop_outer);
        ipr.SetShakeTimes(_ipr_loop_inner);
        ipr.SetTol2D(_ipr_tol_2d);
        ipr.SetTol3D(_ipr_tol_3d);
        ipr.SetShakeWidth(_shake_shift);

        std::vector<T> obj_list;
        ipr.RunIPR(obj_list, _ipr_is_reduced, 1);


        // Link each _active_short_track with a obj candidate
        int n_short_tr = _active_short_track.size();
        std::vector<int> is_erase(n_short_tr,0);
        n_obj = obj_list.size();
        std::vector<int> is_obj_used(n_obj, 0);
        if (n_obj)
        {
            t_start = clock();

            #pragma omp parallel
            {
                #pragma omp for
                for (int i = 0; i < n_short_tr; i ++)
                {
                    MakeShortLinkResidual(nextframe, obj_list, _active_short_track[i], 5, is_erase[i], is_obj_used);
                }
            }

            // Delete labeled tr
            //  from back to front
            for (int i = n_short_tr-1; i > -1; i --)
            {
                if (is_erase[i])
                {
                    _active_short_track.erase(_active_short_track.begin()+i);
                }
            }

            // Delete labeled candidates
            //  from back to front
            for (int i = n_obj-1; i > -1; i --)
            {
                if (is_obj_used[i])
                {
                    obj_list.erase(obj_list.begin() + i);
                }
            }

            t_end = clock();
            std::cout << "Link short track time: " 
                      << (t_end - t_start)/CLOCKS_PER_SEC << "s"
                      << std::endl;
            
            // Move all activeShortTracks longer than 3 particles to activeLongTracks
            //  from back to front
            n_short_tr = _active_short_track.size();
            for (int i = n_short_tr-1; i > -1; i --)
            {
                if (_active_short_track[i].Length() > 3)
                {
                    _active_long_track.push_back(_active_short_track[i]);
                    _active_short_track.erase(_active_short_track.begin()+i);

                    _s_as4 ++;
                    _a_al ++;
                }
            }

            // Add all the untracked candidates to a new short track
            n_obj = obj_list.size();
            for (int i = 0; i < n_obj; i ++)
            {
                _active_short_track.push_back(Track<T> (obj_list[i], nextframe));
            }

            std::cout << "Done." << std::endl;
        }
        

        // Prune / Arrange the tracks
        t_start = clock();

        int n_long_tr = _active_long_track.size();
        std::cout << "debug:n_long_tr=" << n_long_tr << std::endl;
        is_erase = std::vector<int> (n_long_tr, 0);
        #pragma omp parallel
        {
            #pragma omp for
            for (int i = 0; i < n_long_tr; i ++)
            {
                if (!CheckLinearFit(_active_long_track[i]))
                {
                    is_erase[i] = 1;
                }
            }
        }

        // Delete unstatisfied long tr
        //  from back to front
        for (int i = n_long_tr-1; i > -1; i --)
        {
            if (is_erase[i])
            {
                if (_active_long_track[i].Length() > 6)
                {
                    _inactive_long_track.push_back(_active_long_track[i]);
                        
                    _a_il ++;
                }
                else
                {
                    _a_is ++;
                }

                _s_al ++;
                _active_long_track.erase(_active_long_track.begin()+i);
            }
        }
        

        // Print info
        t_end = clock();

        std::cout << "Pruning time: " << (t_start-t_end)/CLOCKS_PER_SEC << "s" << std::endl;

        std::cout << "\t\tNo. of active Short tracks: " << c1 << " + " << _a_as << " - (" << _s_as1 << " + " << _s_as2 << " + " << _s_as3 << " + " << _s_as4 << ") = " << _active_short_track.size() << std::endl;

        std::cout << "\t\tNo. of active Long tracks: " << c2 << " + " << _a_al << " - " << _s_al << " = " << _active_long_track.size() << std::endl;
        
        std::cout << "\t\tNo. of exited tracks: " << _exit_track.size() << std::endl;

        std::cout << "\t\tNo. of inactive Long tracks:	" << c4 << " + " << _a_il << " = " << _inactive_long_track.size() << std::endl;


        // Save data
        //  Save the inacitve long tracks for every 500 frames 
        //  and empty it to avoid endless expansion of this variable
        if (nextframe % 10 == 0)
        {
            // save the inactive long tracks
            std::string s = std::to_string(nextframe);
            std::string X5 = "InactiveLongTrack" + s + ".csv";
            SaveTrack(_inactive_long_track, address + X5);
            // empty inactiveLongTracks
            if (nextframe != endframe) 
            {
                _inactive_long_track.erase(_inactive_long_track.begin(), _inactive_long_track.end());
            }

            // save the inactive long tracks
            std::string X6 = "ExitTrack" + s + ".csv";
            SaveTrack(_exit_track, address + X6);
            // empty inactiveLongTracks
            if (nextframe != endframe) 
            {
                _exit_track.erase(_exit_track.begin(), _exit_track.end());
            }
        }
    }


    // Save all the data
    std::string s = std::to_string(endframe-1);
    std::string file;

    // save the active long tracks
    file = "ActiveLongTrack" + s + ".csv";
    SaveTrack(_active_long_track, address + file);

    // save the active short tracks
    file = "ActiveShortTrack" + s + ".csv";
    SaveTrack(_active_short_track, address + file);

    //  save the inactive tracks
    file = "InactiveTrack" + s + ".csv";
    SaveTrack(_inactive_track, address + file);

    //  save the exit tracks
    file = "ExitTrack" + s + ".csv";
    SaveTrack(_exit_track, address + file);

    // save the inactive long tracks
    file = "InactiveLongTrack" + s + ".csv";
    SaveTrack(_inactive_long_track, address + file);

}


template<class T>
void STB<T>::Run ()
{
    // Create folder    
    // create track folder
    CreateFolder(_img_folder + "Track/");
    
    // create initial track folder
    CreateFolder(_img_folder + "Track/InitialTrack/");

    // create converged track folder
    CreateFolder(_img_folder + "Track/ConvergedTrack/");
    

    // Start STB
    InitialPhase();

    clock_t t_start, t_end;
    t_start = clock();

    ConvergencePhase();

    t_end = clock();

    std::cout << "\n\nTotal time taken by STB convergence phase: " << (t_end-t_start)/CLOCKS_PER_SEC << "s" << std::endl;
}


template<class T>
void STB<T>::MakeLink(int nextframe, const Matrix<double>& vel_curr, double radius, Track<T>& track, bool& active)
{
    Matrix<double> pt_last(3,1);
    Matrix<double> pt_estimate(3,1);

    pt_last = track.Last().GetCenterPos();
    pt_estimate = pt_last + vel_curr;

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
        if (obj_list[i].IsTrack())
        {
            continue;
        }

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
void STB<T>::NearestNeighbor(std::vector<T>& obj_list, double radius, const Matrix<double>& pt_estimate, int& obj_id, std::vector<int>& candidate_used)
{
    obj_id = -1;
    double min2 = radius * radius;
    double dis2 = 0;

    Matrix<double> pt(3,1);
    for (int i = 0; i < obj_list.size(); i ++)
    {
        if (candidate_used[i])
        {
            continue;
        }
        
        pt = obj_list[i].GetCenterPos();
        // dis2 = pt_estimate.DistSqr(pt);
        dis2 = pt.DistSqr(pt_estimate);

        if (dis2 < min2)
        {
            obj_id = i;
            min2 = dis2;
        }
    }

    if (obj_id != _UNLINKED)
    {
        candidate_used[obj_id] = 1;
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


template<class T>
void STB<T>::Prediction(int frame, std::vector<Matrix<double>>& est_pos)
{
    // Use polynomial fit / Wiener filter to predict the particle position at nextFrame
    std::vector<std::vector<double>> pred_coeff(3);

    std::vector<std::string> direction = {"X", "Y", "Z"};
    Matrix<double> est(3,1);

    for (int k = _active_long_track.size()-1; k > -1; k --)
    {
        int size = _active_long_track[k].Length();
        for (int i = 0; i < 3; i ++)
        {
            if (size < 4)
            {
                est(i,0) = LMSWienerPred(_active_long_track[k], direction[i], 3);
            }
            else if (size < 6)
            {
                est(i,0) = LMSWienerPred(_active_long_track[k], direction[i], size-1);
            }
            else 
            {
                est(i,0) = LMSWienerPred(_active_long_track[k], direction[i], 5);
            }
        }

        // Check the boundary
        // if the estimate particle is inside the boundary
        if (_xyz_limit.Check(est(0,0), est(1,0), est(2,0)))
        {
            est_pos.push_back(est);  
        }
        else
        {
            // the track should be put into exit tracks since it is outside the boundary
            _exit_track.push_back(_active_long_track[k]);
            _active_long_track.erase(_active_long_track.begin()+k);
            _s_al ++;
        }
    }
}


template<class T>
double STB<T>::LMSWienerPred(Track<T>& track, std::string direction, int order)
{
    int size = track.Length();
    std::vector<double> series (order+1, 0);

    // Get the series to be predicted
    int axis = -1;
    if (direction == "X" || direction == "x")
    {
        axis = 0;
    }
    else if (direction == "Y" || direction == "y")
    {
        axis = 1;
    }
    else if (direction == "Z" || direction == "z")
    {
        axis = 2;
    }
    else
    {
        std::cout << "STB::LMSWienerPred: no such direction: " 
                    << direction
                    << std::endl;
        throw error_type;
    }   
    for (int i = 0; i < order+1; i ++)
    {
        series[i] = track.GetPos(size-order-1 + i)(axis,0);
    }

    // Wiener filter does badly near zero
    // Make a shift to avoid zero-plane prediction
    bool shift_label = false;
    double shift = 10;
    if (std::fabs(series[order]) < 1)
    {
        shift_label = true;
        for (int i = 0; i < order+1; i ++)
        {
            series[i] = series[i] + shift;
        }
    }

    std::vector<double> filter_param (order, 0);
    
    // calculate the step
    double sum = 0;
    for (int i = 0; i < order; i ++)
    {
        sum += series[i] * series[i];
    }
    double step = 1 / sum;

    double prediction = 0;
    for (int i = 0; i < order; i ++)
    {
        prediction += (filter_param[i] * series[i]);
    }

    double error = series[order] - prediction;

    while (std::fabs(error) > SMALLNUMBER)
    {
        for (int i = 0; i < order; i ++)
        {
            filter_param[i] += (step * series[i] * error);
        }

        // Calculate the prediction using the new filter param
        prediction = 0;
        for (int i = 0; i < order; i ++)
        {
            prediction += filter_param[i] * series[i];
        }
        error = series[order] - prediction;
    }
    prediction = 0;
    for (int i = 0; i < order; i ++)
    {
        prediction += (filter_param[i] * series[i+1]);
    }
    if (shift_label)
    {
        prediction -= shift;
    }

    return prediction;
}


// Link short tracks with obj candidates in residual images
template<class T>
void STB<T>::MakeShortLinkResidual(int nextframe, std::vector<T>& obj_list, Track<T>& tr, int n_iter, int& is_erase, std::vector<int>& candidate_used)
{
    // iteratively trying to find a link for the short track from particle candidates
    Matrix<double> vel(3,1,0);
    Matrix<double> est(3,1);
    int n_long_tr = _active_long_track.size();
    double w;
    double beta2;
    int obj_id = _UNLINKED;

    for (int i = 0; i < n_iter; i ++)
    {
        double rsqr = std::pow(std::pow(1.1, i) * 3 * _space_avg, 2);
        double shift = std::pow(1.1, i) * _shift_betframe_max;
        std::vector<double> weight;
        double tot_weight = 0;
        std::vector<Matrix<double>> disp;
        vel *= 0;
        beta2 = 1 / std::pow(std::pow(1.1, i) * 0.1 * _space_avg, 2); // typical length scale
        
        // calculate the predictive vel. field as an avg. of particle vel from neighbouring tracks
        // identify the neighbour tr (using 3*avg interpt dist.) and get their pt vel
        for (int j = 0; j < n_long_tr; j ++)
        {
            double dsqr = tr.Last().GetCenterPos().DistSqr(_active_long_track[j].Penultimate().GetCenterPos());

            if (dsqr < rsqr && dsqr > 0)
            {
                w = 1/(1 + dsqr*beta2);
                tot_weight = tot_weight + w;
                weight.push_back(w);
                disp.push_back(_active_long_track[j].Last().GetCenterPos()-_active_long_track[j].Penultimate().GetCenterPos());
            }
        }

        int n_neighbour = weight.size();
        if (n_neighbour > 0)
        {
            for (int j = 0; j < n_neighbour; j ++)
            {
                w = weight[j] / tot_weight;              
                vel(0,0) += w * disp[j](0,0);
                vel(1,0) += w * disp[j](1,0);
                vel(2,0) += w * disp[j](2,0);
            }

            est = tr.Last().GetCenterPos() + vel;
            NearestNeighbor(obj_list, _r_search_init, est, obj_id, candidate_used);
        }
        else
        {
            // if no neighbouring tracks are identified
            est = tr.Last().GetCenterPos();
            NearestNeighbor(obj_list, shift, est, obj_id, candidate_used);
        }

        // If linked with a candidate
        //  add the candidate to the track
        if (obj_id != _UNLINKED)
        {
            tr.AddNext(obj_list[obj_id], nextframe);
            break;
        }
    }

    // If no link is found for the short tr
    if (obj_id == _UNLINKED)
    {
        int length = tr.Length();
        if (length > 3)
        {
            // if the track has at least 4 particles
            // _inactive_track.push_back(tr);
            _a_is ++;
        }
        else if (length == 1)
        {
            _s_as1 ++;
        }
        else if (length == 2)
        {
            _s_as2++;
        }
        else if (length == 3)
        {
            _s_as3++;
        }

        is_erase = 1;
    }

}


template<class T>
bool STB<T>::CheckLinearFit(Track<T>& tr)
{
    int len = tr.Length();

    std::vector<Matrix<double>> pt_list (4, Matrix<double> (3,1,0));
    for (int i = 0; i < 4; i ++)
    {
        pt_list[i] = tr.GetPos(len - 4 + i); 
    }

    double y0,y1,y2,y3;
    Matrix<double> coeff(3,2);
    Matrix<double> est(3,1,0);
    for (int i = 0; i < 3; i ++)
    {
        y0 = pt_list[0](i,0);
        y1 = pt_list[1](i,0);
        y2 = pt_list[2](i,0);
        y3 = pt_list[3](i,0);

        coeff(i,0) = - 0.3*y0 - 0.1*y1 + 0.1*y2 + 0.3*y3;
        coeff(i,1) =   0.7*y0 + 0.4*y1 + 0.1*y2 - 0.2*y3;
        est(i,0) = coeff(i,0) * 3 + coeff(i,1);
    }

    double dist[4];
    dist[3] = est.Dist(pt_list[3]);

    std::cout << dist[3] << std::endl;

    double err = 0.05; // 0.05 mm
    if (dist[3] > err) 
    {
        return false;
    }
    else if (dist[3] > err * 0.6)
    {
        // For those tracks that are at the edge of being deleted, check the whole track to see whether its error is large also
        double tot_dist = dist[3];
        for (int i = 0; i < 3; i ++)
        {
            est(0,0) = coeff(0,0) * i + coeff(0,1);
            est(1,0) = coeff(1,0) * i + coeff(1,1);
            est(2,0) = coeff(2,0) * i + coeff(2,1);

            dist[i] = est.Dist(pt_list[i]);
            if (dist[i] > err * 0.6)
            {
                return false;
            }
            tot_dist += dist[i];
        }
        if (tot_dist/4 > err * 0.6)
        {
            return false;
        }
    }

    return true;
}




// DATA IO //
template<class T>
void STB<T>::CreateFolder(std::string folder)
{
    int ret = mkdir(folder.c_str());
    if (ret && errno==EEXIST)
    {
        std::cout << "DIR: " << folder << " already exists!" << std::endl;
    }
    else if (ret)
    {
        std::cout << "Create DIR error: " << ret << strerror(errno) << std::endl;
        throw;
    }
    else
    {
        std::cout << "DIR is successfully created!" << std::endl;
    }
}

template<class T>
void STB<T>::SaveTrack(std::deque<Track<T>> tracks, std::string address)
{
    int n_tr = tracks.size();

    int tot_len = 0; // total length of all the tracks
    for (int i = 0; i < n_tr; i ++)
    {
        tot_len += tracks[i].Length();
    }

    std::ofstream outfile(address, std::ios::out);

    Matrix<double> pt(3,1);
    for (int i = 0; i < n_tr; i ++)
    {
        int len = tracks[i].Length();
        int t_start = tracks[i].GetTime(0);

        for (int j = 0; j < len; j ++)
        {
            pt = tracks[i].GetPos(j);
            outfile << i << ","; // track ID
            outfile << (t_start+j) << ","; // frame ID
            outfile << pt(0,0) << "," << pt(1,0) << "," << pt(2,0) << "\n"; // 3D pos
        }
    }

    outfile.close();
}


// PENDING
template<class T> 
void STB<T>::LoadTrack(std::string file, TrackType trackType)
{

}

#endif