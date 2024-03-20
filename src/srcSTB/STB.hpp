#ifndef STB_HPP
#define STB_HPP

#include "STB.h"

template<class T3D>
STB<T3D>::STB(int frame_start, int frame_end, float fps, double vx_to_mm, int n_thread, std::string const& output_folder, CamList const& cam_list, AxisLimit const& axis_limit,  std::string const& file)
    : _first(frame_start), _last(frame_end), _fps(fps), _vx_to_mm(vx_to_mm), _n_thread(n_thread), _output_folder(output_folder), _cam_list(cam_list), _n_cam_all(cam_list.cam_list.size()), _axis_limit(axis_limit)
{
    // Create output folder
    createFolder(output_folder);
    createFolder(output_folder + "InitialTrack/");
    createFolder(output_folder + "ConvergeTrack/");

    // Read param
    std::ifstream stb_config(file, std::ios::in);

    if (!stb_config.is_open())
    {
        std::cerr << "STB<T3D>::STB error at line" << __LINE__ << ":\n"
                  << "Cannot open file " << file << std::endl;
        throw error_io;
    }

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
        }
        else if (commentpos == 0)
        {
            continue;
        }

        parsed << line << '\t';
    }
    stb_config.close();
    

    // Load Tracking parameters //
    parsed >> _ipr_flag; // 1 for using ipr in initial phase (0 for use .csv files)
    parsed >> _r_objSearch; _r_objSearch *= _vx_to_mm; // mm
    parsed >> _n_initPhase; // number of frames for initial phase
    parsed >> _r_trackSearch; _r_trackSearch *= _vx_to_mm;

    // Load Shake parameters //
    parsed >> _shake_width; _shake_width *= _vx_to_mm;

    // Load Predict Field parameters //
    _pf_param.limit = _axis_limit;
    parsed >> _pf_param.nx >> _pf_param.ny >> _pf_param.nz;
    parsed >> _pf_param.r; _pf_param.r *= _vx_to_mm;

    // Load IPR parameters //
    _ipr_param.n_thread = _n_thread;

    int option;
    parsed >> option;
    if (option == 1)
    {
        _ipr_param.tri_only = true;
    }
    else if (option == 2)
    {
        _ipr_only = true;
    }

    parsed >> _ipr_param.n_loop_ipr;
    parsed >> _ipr_param.n_loop_shake;
    parsed >> _ipr_param.ghost_threshold;
    parsed >> _ipr_param.tol_2d; // [px]
    parsed >> _ipr_param.tol_3d; _ipr_param.tol_3d *= _vx_to_mm;
    parsed >> _n_reduced;
    parsed >> _ipr_param.n_loop_ipr_reduced;

    // std::cout << "IPR parameters: " << std::endl;
    // std::cout << "\tNumber of IPR loop: " << _ipr_param.n_loop_ipr << std::endl;
    // std::cout << "\tNumber of IPR loop for shake: " << _ipr_param.n_loop_shake << std::endl;
    // std::cout << "\tGhost threshold: " << _ipr_param.ghost_threshold << std::endl;
    // std::cout << "\tTolerance for 2D matching: " << _ipr_param.tol_2d << " [px]" << std::endl;
    // std::cout << "\tTolerance for 3D matching: " << _ipr_param.tol_3d << " [mm]" << std::endl;
    // std::cout << "\tNumber of reduced camera: " << _n_reduced << std::endl;
    // std::cout << "\tNumber of IPR loop for reduced camera: " << _ipr_param.n_loop_ipr_reduced << std::endl;


    // Load Object Info //
    loadObjParam(parsed);

    std::cout << std::endl;
}


template<class T3D>
void STB<T3D>::processFrame (int frame_id, std::vector<Image>& img_list)
{
    if (frame_id < _first || frame_id > _last)
    {
        std::cout << "Frame " << frame_id << " is out of range!" << std::endl;
        return;
    }

    // Process STB
    clock_t t_start, t_end;
    t_start = clock();

    if (frame_id-_first < _n_initPhase)
    {
        // Initial phase
        std::cout << "Initial phase at frame " << frame_id << std::endl;

        runInitPhase(frame_id, img_list);
    }
    else
    {
        // Convergence phase
        std::cout << "Convergence phase at frame " << frame_id << std::endl;

        runConvPhase(frame_id, img_list);
    }

    t_end = clock();
    std::cout << "Total time for frame " << frame_id << ": " << (double) (t_end - t_start)/CLOCKS_PER_SEC << std::endl;
    std::cout << std::endl;

    if (frame_id == _last)
    {
        std::cout << "Finish STB!" << std::endl;
        std::cout << "Frame start: " << _first << "; Frame end: " << _last << std::endl;
        std::cout << "Number of initial phase frames: " << _n_initPhase << std::endl;

        std::cout << "Output folder: " << _output_folder << std::endl;

        std::cout << "Number of active short tracks: " << _short_track_active.size() << std::endl;
        std::cout << "Number of active long tracks: " << _long_track_active.size() << std::endl;
        std::cout << "Number of inactive long tracks: " << _long_track_inactive.size() << std::endl;
        std::cout << "Number of exited tracks: " << _exit_track.size() << std::endl;

        saveTracksAll(_output_folder, frame_id);
    }
}


template<class T3D>
void STB<T3D>::createFolder (std::string const& folder)
{
    // // Check if src folder exists
    // if (!fs::is_directory(folder) || !fs::exists(folder)) 
    // { 
    //     fs::create_directories(folder); // create src folder
    // }

    fs::create_directories(folder);
}


template<class T3D>
void STB<T3D>::loadObjParam (std::stringstream& config)
{
    if (typeid(T3D) == typeid(Tracer3D))
    {
        _obj_param.resize(3);
        for (int i = 0; i < 3; i ++)
        {
            config >> _obj_param[i];
        }

        // initialize OTF
        _otf.loadParam(_n_cam_all, 2, 2, 2, _axis_limit);

        // std::cout << "Object parameters: " << std::endl;
        // std::cout << "\tObject radius: " << _obj_param[0] << " [px]" << std::endl;
        // std::cout << "\tObject max intensity: " << _obj_param[1] << std::endl;
        // std::cout << "\tObject min intensity: " << _obj_param[2] << std::endl;
        // std::cout << std::endl;
    }
    else
    {
        std::cerr << "STB<T3D>::loadObjParam error!" << std::endl;
        throw error_type;
    }
}


template<class T3D>
void STB<T3D>::runInitPhase (int frame, std::vector<Image>& img_list)
{
    // Initialize IPR
    IPR ipr(_cam_list, img_list, _ipr_param);

    if (_ipr_flag)
    {
        // IPR
        std::vector<T3D> obj3d_list;
        _ipr_matched.push_back(obj3d_list);
        
        ipr.runIPR(_ipr_matched.back(), _obj_param, _otf, _n_reduced);
    }
    else
    {
        // TODO: Load .csv files
    }

    if (frame-_first == _n_initPhase-1)
    {
        if (!_ipr_only && !_ipr_param.tri_only)
        {
            std::cout << std::endl;

            Pt3D vel_curr;
            for (int i = 0; i < _n_initPhase-1; i ++)
            {   
                int currframe = _first + i;
                int nextframe = currframe + 1;

                std::cout << "STB initial phase tracking b/w frames: "
                          << currframe << " & " << nextframe << "; ";
                
                // Predict Field
                PredField pf(_pf_param, _ipr_matched[i], _ipr_matched[i+1]);
                pf.saveDispField(_output_folder + "PredField_" + std::to_string(currframe) + "_" + std::to_string(nextframe) + ".csv");

                // Link particles b/w two frames
                clock_t t_start, t_end;
                t_start = clock();

                // Extend all the tracks that are active in current frame
                std::vector<int> is_tracked_next(_ipr_matched[i+1].size(), 0);
                for (typename std::deque<Track<T3D>>::iterator track = _short_track_active.begin(); track != _short_track_active.end();)
                {
                    pf.getDisp(vel_curr, track->_obj3d_list.back()._pt_center);
                    LinkStatus status = makeLink(*track, nextframe, vel_curr, _r_objSearch);

                    if (! status.active)
                    {
                        track = _short_track_active.erase(track);
                    }
                    else
                    {
                        ++ track;
                        is_tracked_next[status.obj3d_id] = 1;
                    }
                }

                // Update _is_tracked status for the next frame
                for (int j = 0; j < _ipr_matched[i+1].size(); j ++)
                {
                    if (is_tracked_next[j])
                    {
                        _ipr_matched[i+1][j]._is_tracked = true;
                    }
                }

                // Start a track for all particles left untracked in current frame
                startTrack(currframe, pf);

                t_end = clock();
                std::cout << "; total link particles time: " 
                        << (double) (t_end - t_start)/CLOCKS_PER_SEC
                        << std::endl;
            }

            // Move all tracks longer than 3 from _short_track_active to _long_track_active
            for (typename std::deque<Track<T3D>>::iterator track = _short_track_active.begin(); track != _short_track_active.end();)
            {
                if (track->_obj3d_list.size() > 3)
                {
                    _long_track_active.push_back(*track);
                    track = _short_track_active.erase(track);
                }
                else
                {
                    ++ track;
                }
            }

            // Add the untracked particles at the last initial phase frame (_first+_n_initPhase-1) to _short_track_active
            int m = _n_initPhase - 1;
            for (int i = 0; i < _ipr_matched[m].size(); i ++)
            {
                if (!_ipr_matched[m][i]._is_tracked)
                {
                    _short_track_active.push_back(Track<T3D> (_ipr_matched[m][i], _first + m));
                }
            }

            std::cout << "Done with initial phase!" << std::endl;
            std::cout << "\tNo. of active short tracks: " << _short_track_active.size() 
                    << "; No. of active long tracks: " << _long_track_active.size()
                    << "; No. of exited tracks: " << _exit_track.size()
                    << std::endl;
        }


        // Save all data
        std::string address = _output_folder + "InitialTrack/";

        // ipr data
        if (_ipr_flag)
        {
            for (int i = 0; i < _n_initPhase; i ++)
            {
                ipr.saveObjInfo(address + "IPR_" + std::to_string(_first+i) + ".csv", _ipr_matched[i]);
            }
        }
        
        // save all tracks
        saveTracksAll(address, frame);
    }

}


template<class T3D>
void STB<T3D>::runConvPhase (int frame, std::vector<Image>& img_list)
{
    
}


// TODO: I donot update _is_tracked on the fly (to avoid conflict during parallelization), not sure whether this will affect the result
// Note: I still update the _is_track status after makeLink, but set all possible links to be true
template<class T3D>
LinkStatus STB<T3D>::makeLink(Track<T3D>& track, int nextframe, Pt3D const& vel_curr, double radius)
{
    LinkStatus status;
    status.active = false;
    
    int m = nextframe-_first;

    int obj_id = findNN(_ipr_matched[m], track._obj3d_list.back()._pt_center+vel_curr, radius);

    if (obj_id != UNLINKED)
    {
        track.addNext(_ipr_matched[m][obj_id], nextframe);
        status.active = true;
        status.obj3d_id = obj_id;
    }

    return status;
}


template<class T3D>
int STB<T3D>::findNN(std::vector<T3D> const& obj3d_list, Pt3D const& pt3d_est, double radius)
{
    double min_dist = radius;
    int obj_id = UNLINKED;

    double dist;
    for (int i = 0; i < obj3d_list.size(); i ++)
    {
        if (obj3d_list[i]._is_tracked)
        {
            continue;
        }

        dist = myMATH::dist(obj3d_list[i]._pt_center, pt3d_est);
        if (dist < min_dist)
        {
            min_dist = dist;
            obj_id = i;
        }
    }

    return obj_id;
}   


template<class T3D>
void STB<T3D>::startTrack (int frame, PredField& pf)
{
    int m = frame - _first;
    int n_obj3d = _ipr_matched[m].size();
    int n_obj3d_next = _ipr_matched[m+1].size();

    std::vector<int> is_tracked_next(n_obj3d_next, 0);
    for (int i = 0; i < n_obj3d; i ++)
    {
        if (!_ipr_matched[m][i]._is_tracked)
        {
            _ipr_matched[m][i]._is_tracked = true;

            // Start a track for the untracked particle            
            Track<T3D> init_tr(_ipr_matched[m][i], frame);

            Pt3D vel_curr;
            pf.getDisp(vel_curr, _ipr_matched[m][i]._pt_center);
            LinkStatus status = makeLink(init_tr, frame+1, vel_curr, _r_objSearch);

            if (status.active)
            {
                _short_track_active.push_back(init_tr);
                is_tracked_next[status.obj3d_id] = 1;
            }
        }
    }

    // Update _is_tracked status for the next frame
    for (int i = 0; i < n_obj3d_next; i ++)
    {
        if (is_tracked_next[i])
        {
            _ipr_matched[m+1][i]._is_tracked = true;
        }
    }
}


template<class T3D>
void STB<T3D>::saveTracks (std::string const& file, std::deque<Track<T3D>>& tracks)
{
    std::ofstream output(file, std::ios::out);

    if (!output.is_open())
    {
        std::cerr << "STB<T3D>::saveTracks error at line" << __LINE__ << ":\n"
                  << "Cannot open file " << file << std::endl;
        throw error_io;
    }

    output << "track_id,t,world_x,world_y,world_z,error";
    for (int i = 0; i < _n_cam_all; i ++)
    {
        output << ",cam" << i << "_x(col),cam" << i << "_y(row)";
    }
    output << "\n";

    for (int i = 0; i < tracks.size(); i ++)
    {
        tracks[i].saveTrack(output, i, _fps, _n_cam_all);
    }

    output.close();
}


template<class T3D>
void STB<T3D>::saveTracksAll(std::string const& folder, int frame)
{
    std::string s = std::to_string(frame);
    std::string file;

    file = "LongTrackActive_" + s + ".csv";
    saveTracks(folder + file, _long_track_active);

    file = "ShortTrackActive_" + s + ".csv";
    saveTracks(folder + file, _short_track_active);

    file = "ExitTrack_" + s + ".csv";
    saveTracks(folder + file, _exit_track);

    file = "LongTrackInactive_" + s + ".csv";
    saveTracks(folder + file, _long_track_inactive);
}

#endif