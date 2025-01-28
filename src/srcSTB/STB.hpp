#ifndef STB_HPP
#define STB_HPP

#include "STB.h"

// Auxiliary functions //

// Calibrate OTF parameters
// Gaussian intensity: 
//  projection center: (xc,yc), x=col, y=row
//  dx =  (x-xc)*cos(alpha) + (y-yc)*sin(alpha)
//  dy = -(x-xc)*sin(alpha) + (y-yc)*sin(alpha)
//  I(x,y) = a * exp( - b*dx^2 - c*dx^2 )
//  - ln(I) + ln(a) = b*dx^2 + c*dy^2
//  a = mean(I_c) + 2*std(I_c)
// Assume a,b,c,alpha are the same in the view volume, alpha=0
//  nx,ny,nz = 2,2,2
// otf_param: need to be initialized before calling this function
void calibOTFParam (OTFParam& otf_param, int cam_id, double radius, std::vector<std::vector<Tracer2D>> const& tr2d_list, std::vector<Image> const& imgList)
{
    // check if the camera is valid
    if (cam_id < 0 || cam_id >= otf_param.a.getDimRow())
    {
        std::cerr << "calibOTFParam error at line " << __LINE__ << ":\n" << "Invalid camera id: " << cam_id << "; total number of cameras: " << otf_param.a.getDimRow() << std::endl;
        throw error_size;
    }

    // check if the number of tracers is the same as the number of images
    if (tr2d_list.size() != imgList.size())
    {
        std::cerr << "calibOTFParam error at line " << __LINE__ << ":\n" << "The number of tracer lists is not the same as the number of images!" << "Number of tracer list: " << tr2d_list.size() << ". Number of img: " << imgList.size() << std::endl;
        throw error_size;
    }

    std::vector<double> I_list;
    std::vector<double> a_list;
    std::vector<std::vector<double>> coeff_list;
    std::vector<double> coeff(2,0);
    double xc, yc;
    int n_row, n_col;
    int row_min, row_max, col_min, col_max;
    for (int i = 0; i < imgList.size(); i ++)
    {
        n_row = imgList[i].getDimRow();
        n_col = imgList[i].getDimCol();

        for (int j = 0; j < tr2d_list[i].size(); j ++)
        {
            xc = tr2d_list[i][j]._pt_center[0];
            yc = tr2d_list[i][j]._pt_center[1];

            a_list.push_back(imgList[i](std::round(yc), std::round(xc)));

            row_min = std::max(0, int(std::floor(yc - radius)));
            row_max = std::min(n_row, int(std::ceil(yc + radius + 1)));
            col_min = std::max(0, int(std::floor(xc - radius)));
            col_max = std::min(n_col, int(std::ceil(xc + radius + 1)));

            for (int row = row_min; row < row_max; row ++)
            {
                for (int col = col_min; col < col_max; col ++)
                {
                    I_list.push_back(imgList[i](row, col));
                    coeff[0] = (std::pow(col-xc, 2));
                    coeff[1] = (std::pow(row-yc, 2));
                    coeff_list.push_back(coeff);
                }
            }
        }
    }

    // a = mean(I_c) + 2.*std(I_c)
    double a;
    double a_mean = 0;
    double a_std = 0;
    double a_max = *std::max_element(a_list.begin(), a_list.end());
    int n_a = a_list.size();
    for (int i = 0; i < n_a; i ++)
    {
        a_mean += a_list[i];
    }
    a_mean /= n_a;
    for (int i = 0; i < n_a; i ++)
    {
        a_std += std::pow(a_list[i] - a_mean, 2);
    }
    a_std = std::sqrt(a_std / n_a);
    a = std::min(a_mean + 2 * a_std, a_max);  
    std::cout << "\ta_mean = " << a_mean << "; a_std = " << a_std << "; a_max = " << a_max << "; a = " << a << std::endl;

    // estimate the coefficients
    int n_coeff = coeff_list.size();
    double logI;
    Matrix<double> coeff_mat(n_coeff, 2, 0);
    Matrix<double> logI_mat(n_coeff, 1, 0);
    for (int i = 0; i < n_coeff; i ++)
    {
        for (int j = 0; j < 2; j ++)
        {
            coeff_mat(i, j) = coeff_list[i][j];
        }

        logI = I_list[i] < LOGSMALLNUMBER ? std::log(LOGSMALLNUMBER) : std::log(I_list[i]);
        logI_mat(i, 0) = - logI + std::log(a);
    }
    Matrix<double> coeff_est = myMATH::inverse(coeff_mat.transpose() * coeff_mat) * coeff_mat.transpose() * logI_mat;

    for (int i = 0; i < 8; i ++)
    {
        otf_param.a(cam_id, i) = a;
        otf_param.b(cam_id, i) = coeff_est(0, 0);
        otf_param.c(cam_id, i) = coeff_est(1, 0);
        otf_param.alpha(cam_id, i) = 0;
    }
}


// Member functions of STB //

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
    parsed >> _r_predSearch; // radius to find predicted object [px]

    if (_n_initPhase < 2)
    {
        std::cerr << "STB<T3D>::STB error at line" << __LINE__ << ":\n"
                  << "Number of frames for initial phase should be at least 2!" << std::endl;
        throw error_size;
    }

    // load tracking predictor parameters
    parsed >> line; 
    bool is_pred_param = false;
    if (line.rfind("Wiener", 0) == 0)
    {
        _track_pred_param.type = WIENER;
        is_pred_param = true;
    }
    else if (line.rfind("Kalman", 0) == 0)
    {
        _track_pred_param.type = KALMAN;
        // load float param based on comma
        std::stringstream param_ss(line.substr(7));
        double param_i;
        while (param_ss >> param_i)
        {
            param_i *= _vx_to_mm;
            _track_pred_param.param.push_back(param_i);
            if (param_ss.peek() == ',')
            {
                param_ss.ignore();
            }
        }
        is_pred_param = true;
    }

    // Load Shake parameters //
    if (is_pred_param)
    {
        parsed >> _shake_width; _shake_width *= _vx_to_mm;
    }
    else
    {
        _shake_width = std::stod(line); _shake_width *= _vx_to_mm;
    }
    _ipr_param.shake_width = _shake_width;

    // Load Predict Field parameters //
    _pf_param.limit = _axis_limit;
    parsed >> _pf_param.nx >> _pf_param.ny >> _pf_param.nz;
    parsed >> _pf_param.r; _pf_param.r *= _vx_to_mm;
    _pf_param.nBin_x = 11;
    _pf_param.nBin_y = 11;
    _pf_param.nBin_z = 11;

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
    parsed >> _ipr_param.n_obj2d_max; // maximum number of tracers in each camera
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
STB<T3D>::STB(const STB& stb)
    : _ipr_matched(stb._ipr_matched), _short_track_active(stb._short_track_active), _long_track_active(stb._long_track_active), _long_track_inactive(stb._long_track_inactive), _exit_track(stb._exit_track), _first(stb._first), _last(stb._last), _fps(stb._fps), _vx_to_mm(stb._vx_to_mm), _n_thread(stb._n_thread), _output_folder(stb._output_folder), _cam_list(stb._cam_list), _n_cam_all(stb._n_cam_all), _axis_limit(stb._axis_limit), _ipr_flag(stb._ipr_flag), _r_objSearch(stb._r_objSearch), _r_trackSearch(stb._r_trackSearch), _n_initPhase(stb._n_initPhase), _r_predSearch(stb._r_predSearch), _track_pred_param(stb._track_pred_param), _shake_width(stb._shake_width), _pf_param(stb._pf_param), _ipr_only(stb._ipr_only), _ipr_param(stb._ipr_param), _n_reduced(stb._n_reduced), _tol_2d_overlap(stb._tol_2d_overlap), _obj_param(stb._obj_param), _otf(stb._otf), _a_sa(stb._a_sa), _a_la(stb._a_la), _s_sa(stb._s_sa), _s_la(stb._s_la), _a_li(stb._a_li)
{
    // Create output folder
    createFolder(_output_folder);
    createFolder(_output_folder + "InitialTrack/");
    createFolder(_output_folder + "ConvergeTrack/");
}


template<class T3D>
std::vector<double> STB<T3D>::getObjParam() const
{
    return _obj_param;
}


template<class T3D>
void STB<T3D>::calibrateOTF(int cam_id, int n_obj2d_max, double r_otf_calib, std::vector<Image> const& img_list)
{
    if (cam_id < 0 || cam_id >= _n_cam_all)
    {
        std::cerr << "STB<T3D>::calibrateOTF error at line" << __LINE__ << ":\n"
                  << "Camera ID " << cam_id << " is out of range: " << "0 ~ " << _n_cam_all-1 << std::endl;
        throw error_range;
    }

    if (typeid(T3D) == typeid(Tracer3D))
    {
        std::vector<std::vector<Tracer2D>> tr2d_list_all;
        ObjectFinder2D objfinder;

        std::cout << "Camera " << cam_id << std::endl;
        std::cout << "\tNumber of found tracers in each image: ";
        for (int j = 0; j < img_list.size(); j ++)
        {
            std::vector<Tracer2D> tr2d_list;
            objfinder.findObject2D(tr2d_list, img_list[j], _obj_param);

            std::cout << tr2d_list.size();

            // if tr2d_list is too large, randomly select some tracers
            int seed = 123;
            if (tr2d_list.size() > n_obj2d_max)
            {
                std::shuffle(tr2d_list.begin(), tr2d_list.end(), std::default_random_engine(seed));
                tr2d_list.erase(tr2d_list.begin()+n_obj2d_max, tr2d_list.end());

                std::cout << "(" << n_obj2d_max << ")";
            }
            else if (tr2d_list.size() == 0)
            {
                std::cerr << "Quit OTF calibration: No tracer found in camera " << cam_id << std::endl;
                throw error_size;
            }

            std::cout << ",";

            tr2d_list_all.push_back(tr2d_list);
        }

        std::cout << std::endl;

        calibOTFParam(_otf._param, cam_id, r_otf_calib, tr2d_list_all, img_list);

        std::cout << "\t(a,b,c,alpha) = " << _otf._param.a(cam_id,0) << "," << _otf._param.b(cam_id,0) << ","  << _otf._param.c(cam_id,0) << "," << _otf._param.alpha(cam_id,0) << std::endl;
    }

    // Save OTF parameters
    _otf.saveParam(_output_folder + "OTF.txt");
}


template<class T3D>
void STB<T3D>::processFrame (int frame_id, std::vector<Image>& img_list, bool is_update_img)
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
        runInitPhase(frame_id, img_list, is_update_img);
    }
    else
    {
        // Convergence phase
        std::cout << "Convergence phase at frame " << frame_id << std::endl;

        runConvPhase(frame_id, img_list, is_update_img);
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
void STB<T3D>::runInitPhase (int frame, std::vector<Image>& img_list, bool is_update_img)
{
    // Initialize IPR
    IPR ipr(_cam_list, img_list, _ipr_param);
    std::vector<T3D> obj3d_list;
    _ipr_matched.push_back(obj3d_list);

    if (_ipr_flag)
    {
        // IPR
        ipr.runIPR(_ipr_matched.back(), _obj_param, _otf, _n_reduced);

        if (is_update_img)
        {
            img_list = ipr._imgRes_list;
        }
    }
    else
    {
        // TODO: Load .csv files

        // if (is_update_img)
        // {
        //     // only calculate residue images
        //     Shake s(_cam_list, _shake_width, _ipr_param.ghost_threshold, _ipr_param.n_loop_shake, _n_thread);
        //     s.runShake(_ipr_matched.back(), _otf, img_list, true);

        //     img_list = s._imgRes_list;
        // }
    }

    if (frame-_first == _n_initPhase-1)
    {
        if (!_ipr_only && !_ipr_param.tri_only)
        {
            std::cout << std::endl;

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

                // extend all the tracks that are active in current frame
                int n_sa = _short_track_active.size();
                std::vector<int> link_id(n_sa, UNLINKED);
                
                if (_n_thread>0)
                {
                    omp_set_num_threads(_n_thread);
                }
                #pragma omp parallel for
                for (int j = 0; j < n_sa; j ++)
                {
                    Pt3D vel_curr;
                    pf.getDisp(vel_curr, _short_track_active[j]._obj3d_list.back()._pt_center);
                    link_id[j] = makeLink(_short_track_active[j], nextframe, vel_curr, _r_objSearch);
                }

                // update short track and _is_tracked status for the next frame
                for (int j = n_sa-1; j > -1; j --)
                {
                    if (link_id[j] != UNLINKED)
                    {
                        _short_track_active[j].addNext(_ipr_matched[i+1][link_id[j]], nextframe);
                        _ipr_matched[i+1][link_id[j]]._is_tracked = true;
                    }
                    else
                    {
                        _short_track_active.erase(_short_track_active.begin()+j);
                    }
                }


                // Start a track for all particles left untracked in current frame
                startTrack(currframe, pf);


                t_end = clock();
                std::cout << "; total link particles time: " 
                        << (double) (t_end - t_start)/CLOCKS_PER_SEC
                        << std::endl;
            }

            // Move all tracks >= _n_initPhase from _short_track_active to _long_track_active
            for (typename std::deque<Track<T3D>>::iterator track = _short_track_active.begin(); track != _short_track_active.end();)
            {
                if (track->_obj3d_list.size() >= _n_initPhase)
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
                    _short_track_active.push_back(Track<T3D> (_ipr_matched[m][i], _first + m, _track_pred_param));
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
        // _ipr_matched.clear();
        
        // save all tracks
        saveTracksAll(address, frame);
    }

}


template<class T3D>
void STB<T3D>::runConvPhase (int frame, std::vector<Image>& img_list, bool is_update_img)
{
    // Initialize some variables
    int n_sa = _short_track_active.size();
    int n_la = _long_track_active.size();
    int n_li = _long_track_inactive.size();
    _a_sa = 0; _a_la = 0; _s_sa = 0; _s_la = 0; _a_li = 0;


    // Prediction for active long tracks //
    std::cout << "\tPrediction: ";
    clock_t t_start, t_end;
    t_start = clock();

    std::vector<T3D> obj3d_list_pred(n_la);
    std::vector<int> is_inRange(n_la, 1);

    if (_n_thread > 0)
    {
        omp_set_num_threads(_n_thread);
    }
    #pragma omp parallel for
    for (int i = 0; i < n_la; i ++)
    {
        _long_track_active[i].predictNext(obj3d_list_pred[i]);
        is_inRange[i] = _axis_limit.check(
            obj3d_list_pred[i]._pt_center[0], 
            obj3d_list_pred[i]._pt_center[1], 
            obj3d_list_pred[i]._pt_center[2]
        );

        // // Project the predicted particle to all cameras
        // //  check if there is an actual particle in the image
        // if (is_inRange[i])
        // {
        //     is_inRange[i] = findPredObj(obj3d_list_pred[i], img_list);
        // }

        // if less than 2 cameras can see the predicted particle remove it
        if (is_inRange[i])
        {
            is_inRange[i] = checkObj2D(obj3d_list_pred[i]);
        }
    }

    // Remove out-of-range tracks
    for (int i = n_la-1; i >= 0; i --)
    {
        if (!is_inRange[i])
        {
            if (_long_track_active[i]._n_obj3d >= LEN_LONG_TRACK)
            {
                _exit_track.push_back(_long_track_active[i]);
            }
            
            _long_track_active.erase(_long_track_active.begin()+i);
            obj3d_list_pred.erase(obj3d_list_pred.begin()+i);
            _s_la ++;
        }
    }

    // Remove repeated prediction
    int n_pred = obj3d_list_pred.size();
    std::vector<int> is_repeat(n_pred, 0);
    findRepeatObj(
        is_repeat, 
        obj3d_list_pred, 
        2*_ipr_param.tol_3d
    );
    for (int i = n_pred-1; i >= 0; i --)
    {
        if (is_repeat[i])
        {
            if (_long_track_active[i]._n_obj3d >= LEN_LONG_TRACK)
            {
                _long_track_inactive.push_back(_long_track_active[i]);
            }

            _long_track_active.erase(_long_track_active.begin()+i);
            obj3d_list_pred.erase(obj3d_list_pred.begin()+i);
            _s_la ++;
        }
    }

    t_end = clock();
    std::cout << (double) (t_end - t_start)/CLOCKS_PER_SEC << " s. Done!" << std::endl;

    // Shake prediction
    int n_fail_shaking = 0;
    Shake s (
        _cam_list, 
        _shake_width,
        _ipr_param.tol_3d * 2, 
        _ipr_param.ghost_threshold, //0.01, 
        _ipr_param.n_loop_shake, 
        _n_thread
    );

    if (obj3d_list_pred.size() > 0)
    {
        t_start = clock();
 
        s.runShake(obj3d_list_pred, _otf, img_list, false);

        t_end = clock();
        std::cout << " Shake prediction: " << (double) (t_end - t_start)/CLOCKS_PER_SEC << " s. ";


        // Add the corrected particle position in nextFrame to its respective track
        // remove the track if the particle is a ghost
        for (int i = obj3d_list_pred.size()-1; i >= 0; i --)
        {
            if (!s._is_ghost[i])
            {
                _long_track_active[i].addNext(obj3d_list_pred[i], frame);
            }
            else
            {
                n_fail_shaking ++;

                if (_long_track_active[i]._n_obj3d >= LEN_LONG_TRACK)
                {
                    _long_track_inactive.push_back(_long_track_active[i]);
                    _a_li ++;
                }

                _long_track_active.erase(_long_track_active.begin()+i);
                _s_la ++;
            }
        }

        std::cout << "Finish updating active long tracks!" << std::endl;
    }
    else 
    {
        s._imgRes_list = img_list;
    }


    // IPR on residue images //
    IPR ipr(_cam_list, s._imgRes_list, _ipr_param);
    std::vector<T3D> obj3d_list;
    ipr.runIPR(obj3d_list, _obj_param, _otf, _n_reduced);

    // remove the particles that are close to the active long tracks
    removeOverlap(obj3d_list);

    // update img_list
    if (is_update_img)
    {
        Shake s_resImg (
            _cam_list, 
            _shake_width,
            _ipr_param.tol_3d * 2, 
            _ipr_param.ghost_threshold, //0.01, 
            _ipr_param.n_loop_shake, 
            _n_thread
        );
        
        s_resImg.runShake(obj3d_list, _otf, s._imgRes_list, true);
        img_list = s_resImg._imgRes_list;
    }

    // Link each _short_track_active to an obj 
    int n_obj3d = obj3d_list.size();
    if (n_obj3d > 0)
    {
        std::cout << " Linking: ";
        t_start = clock();

        std::vector<int> link_id(n_sa, 0);
        std::vector<int> is_obj_used(n_obj3d, 0);

        if (_n_thread>0)
        {
            omp_set_num_threads(_n_thread);
        }
        #pragma omp parallel for
        for (int i = 0; i < n_sa; i ++)
        {
            link_id[i] = linkShortTrack (_short_track_active[i], obj3d_list, 5);
        }

        // update _short_track_active and _is_tracked status 
        for (int i = n_sa-1; i >= 0; i --)
        {
            if (link_id[i] != UNLINKED)
            {
                _short_track_active[i].addNext(obj3d_list[link_id[i]], frame);
                obj3d_list[link_id[i]]._is_tracked = true;
            }
            else
            {
                _short_track_active.erase(_short_track_active.begin()+i);
                _s_sa ++;
            }
        }

        // move all active short tracks to long tracks if they have >= _n_initPhase particles
        int n_short_track_active = _short_track_active.size();
        for (int i = n_short_track_active-1; i >= 0; i --)
        {
            if (_short_track_active[i]._n_obj3d >= _n_initPhase)
            {
                _long_track_active.push_back(_short_track_active[i]);
                _short_track_active.erase(_short_track_active.begin()+i);
                _a_la ++;
                _s_sa ++;
            }
        }

        // add all the untracked candidates to a new short track
        for (int i = 0; i < n_obj3d; i ++)
        {
            if (!obj3d_list[i]._is_tracked)
            {
                _short_track_active.push_back(Track<T3D> (obj3d_list[i], frame, _track_pred_param));
                _a_sa ++;
            }
        }


        t_end = clock();
        std::cout << (double) (t_end - t_start)/CLOCKS_PER_SEC << " s. Done!" << std::endl;
    }


    // Prune and arrange the tracks //
    int n_fail_lf = 0;
    // if (frame - _first > _n_initPhase-1)
    if (1)
    {
        t_start = clock();

        int n_la_new = _long_track_active.size();
        std::vector<int> is_erase(n_la_new, 0);

        if (_n_thread>0)
        {
            omp_set_num_threads(_n_thread);
        }
        #pragma omp parallel for
        for (int i = 0; i < n_la_new; i ++)
        {
            if (!checkLinearFit(_long_track_active[i]))
            {
                is_erase[i] = 1;
            }
        }

        // remove the tracks that are not linear    
        for (int i = n_la_new-1; i >= 0; i --)
        {
            if (is_erase[i])
            {
                if (_long_track_active[i]._obj3d_list.size() >= LEN_LONG_TRACK)
                {
                    _long_track_inactive.push_back(_long_track_active[i]);
                    _a_li ++;
                }

                _long_track_active.erase(_long_track_active.begin()+i);
                _s_la ++;
                n_fail_lf ++;
            }
        }

        t_end = clock();

        std::cout << " Pruning time: " << (t_start-t_end)/CLOCKS_PER_SEC << "s" << std::endl;
    }


    // Print outputs
    std::cout << "\tNo. of active short tracks: " << n_sa << " + " << _a_sa << " - " << _s_sa << " = " << _short_track_active.size() << std::endl;

    std::cout << "\tNo. of active long tracks: " << n_la << " + " << _a_la << " - " << _s_la << " = " << _long_track_active.size() << std::endl;
    
    std::cout << "\tNo. of exited tracks: " << _exit_track.size() << std::endl;

    std::cout << "\tNo. of inactive Long tracks: " << n_li << " + " << _a_li << " = " << _long_track_inactive.size() << std::endl;

    std::cout << "\tNo. of fail shaking intensity: " << n_fail_shaking << std::endl;
    std::cout << "\tNo. of fail linear fit: " << n_fail_lf << std::endl;

    // save all data every 500 frames
    if (frame%500 == 0)
    {
        saveTracksAll(_output_folder+"ConvergeTrack/", frame);
        _long_track_inactive.clear();
        _exit_track.clear();
    }
}


template<class T3D>
void STB<T3D>::removeOverlap(std::vector<T3D>& obj3d_list)
{
    if (typeid(T3D) == typeid(Tracer3D))
    {
        removeOverlapTracer(obj3d_list);    
    }
}


template<class T3D>
void STB<T3D>::removeOverlapTracer(std::vector<Tracer3D>& obj3d_list)
{
    int n_obj3d = obj3d_list.size();

    if (_n_thread > 0)
    {
        omp_set_num_threads(_n_thread);
    }
    #pragma omp parallel for
    for (int i = 0; i < n_obj3d; i ++)
    {
        obj3d_list[i].projectObject2D(_cam_list.useid_list, _cam_list.cam_list);
    }

    std::vector<int> is_overlap(n_obj3d, 0);

    #pragma omp parallel for
    for (int i = 0; i < n_obj3d; i ++)
    {
        for (int j = 0; j < _long_track_active.size(); j ++)
        {
            int length = _long_track_active[j]._obj3d_list.size();
            int n_overlap = 0;
            for (int cam_id = 0; cam_id < _n_cam_all; cam_id ++)
            {
                double dx = _long_track_active[j]._obj3d_list[length-1]._tr2d_list[cam_id]._pt_center[0] - obj3d_list[i]._tr2d_list[cam_id]._pt_center[0];
                double dy = _long_track_active[j]._obj3d_list[length-1]._tr2d_list[cam_id]._pt_center[1] - obj3d_list[i]._tr2d_list[cam_id]._pt_center[1];
                bool judge = std::fabs(dx) <= _tol_2d_overlap && std::fabs(dy) <= _tol_2d_overlap;

                if (judge)
                {
                    n_overlap ++;
                }
            }

            if (n_overlap >= _n_cam_all)
            {
                is_overlap[i] = 1;
                break;
            }
        }
    }

    for (int i = n_obj3d-1; i >= 0; i --)
    {
        if (is_overlap[i])
        {
            obj3d_list.erase(obj3d_list.begin()+i);
        }
    }

    // std::cout << "debug:STB::removeOverlapTracer: " << n_obj3d << " -> " << obj3d_list.size() << std::endl;
}


template<class T3D>
void STB<T3D>::findRepeatObj(std::vector<int>& is_repeat, std::vector<T3D> const& obj3d_list, double tol_3d)
{
    int n_obj = obj3d_list.size();
    
    if (is_repeat.size() != n_obj)
    {
        is_repeat.resize(n_obj);
    }
    std::fill(is_repeat.begin(), is_repeat.end(), 0);

    if (_n_thread > 0)
    {
        omp_set_num_threads(_n_thread);
    }

    double repeat_thres_2 = tol_3d * tol_3d;
    for (int i = 0; i < n_obj-1; i ++)
    {
        if (is_repeat[i])
        {
            continue;
        }

        #pragma omp parallel for
        for (int j = i+1; j < n_obj; j ++)
        {
            if (myMATH::dist2(obj3d_list[i]._pt_center, obj3d_list[j]._pt_center) < repeat_thres_2)
            {
                is_repeat[j] = 1;
            }
        }
    }
}


// TODO: I donot update _is_tracked on the fly (to avoid conflict during parallelization), not sure whether this will affect the result
// Note: I still update the _is_track status after makeLink, but set all possible links to be true
template<class T3D>
int STB<T3D>::makeLink(Track<T3D> const& track, int nextframe, Pt3D const& vel_curr, double radius)
{
    int m = nextframe-_first;

    int obj3d_id = findNN(_ipr_matched[m], track._obj3d_list.back()._pt_center+vel_curr, radius);

    return obj3d_id;
}


template<class T3D>
int STB<T3D>::findNN(std::vector<T3D> const& obj3d_list, Pt3D const& pt3d_est, double radius)
{
    double min_dist2 = radius*radius;
    int obj_id = UNLINKED;

    double dist2;
    for (int i = 0; i < obj3d_list.size(); i ++)
    {
        if (obj3d_list[i]._is_tracked)
        {
            continue;
        }

        dist2 = myMATH::dist2(obj3d_list[i]._pt_center, pt3d_est);
        if (dist2 < min_dist2)
        {
            min_dist2 = dist2;
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
    Pt3D vel_curr;

    for (int i = 0; i < n_obj3d; i ++)
    {
        if (!_ipr_matched[m][i]._is_tracked)
        {
            _ipr_matched[m][i]._is_tracked = true;

            // Start a track for the untracked particle      
            Track<T3D> init_tr(_ipr_matched[m][i], frame, _track_pred_param);

            pf.getDisp(vel_curr, _ipr_matched[m][i]._pt_center);
            int obj3d_id = makeLink(init_tr, frame+1, vel_curr, _r_objSearch);

            if (obj3d_id != UNLINKED)
            {
                init_tr.addNext(_ipr_matched[m+1][obj3d_id], frame+1);  
                _short_track_active.push_back(init_tr);
                _ipr_matched[m+1][obj3d_id]._is_tracked = true;
            }
        }
    }
}


template<class T3D>
bool STB<T3D>::findPredObj (T3D& obj3d, std::vector<Image> const& img_list)
{
    Pt2D pt2d;
    PixelRange region;
    int n_row, n_col;

    std::vector<Line3D> sight3D_list;

    for (int i = 0; i < _n_cam_all; i ++)
    {
        pt2d = _cam_list.cam_list[i].project(obj3d._pt_center);
        n_row = _cam_list.cam_list[i].getNRow();
        n_col = _cam_list.cam_list[i].getNCol();

        // create search region
        region.row_min = pt2d[1] - _r_predSearch;
        region.row_max = pt2d[1] + _r_predSearch + 1;
        region.col_min = pt2d[0] - _r_predSearch;
        region.col_max = pt2d[0] + _r_predSearch + 1;

        region.row_min = std::max(0, region.row_min);
        region.row_max = std::min(n_row, region.row_max);
        region.col_min = std::max(0, region.col_min);
        region.col_max = std::min(n_col, region.col_max);

        // find the object in the image
        ObjectFinder2D objfinder;

        if (typeid(T3D) == typeid(Tracer3D))
        {
            obj3d._tr2d_list.clear();
            obj3d._camid_list.clear();

            std::vector<Tracer2D> obj2d_list;
            objfinder.findObject2D(obj2d_list, img_list[i], _obj_param, region);

            if (obj2d_list.size() == 0)
            {
                continue;
            }

            // find nearest object
            std::vector<double> dist2_list(obj2d_list.size(), 0);

            for (int j = 0; j < obj2d_list.size(); j ++)
            {
                dist2_list[j] = myMATH::dist2(obj2d_list[j]._pt_center, pt2d);
            }

            int min_id = std::min_element(dist2_list.begin(), dist2_list.end()) - dist2_list.begin();

            obj3d._tr2d_list.push_back(obj2d_list[min_id]);
            obj3d._camid_list.push_back(i);
            sight3D_list.push_back(_cam_list.cam_list[i].lineOfSight(obj2d_list[min_id]._pt_center));
        }
        else
        {
            std::cerr << "STB<T3D>::findPredObj error at line " << __LINE__ << std::endl;
            std::cerr << "Object type is not supported!" << std::endl;
            throw error_type;
        }
    }


    // the object should be seen by at least 2 cameras
    if (sight3D_list.size() < 2)
    {
        return false;
    }

    // triangulate the object
    myMATH::triangulation(obj3d._pt_center, obj3d._error, sight3D_list);
    if (obj3d._error > _ipr_param.tol_3d)
    {
        return false;
    }


    return true;
}


template<class T3D>
bool STB<T3D>::checkObj2D (T3D& obj3d)
{
    int n_outrange = 0;
    obj3d.projectObject2D(_cam_list.useid_list, _cam_list.cam_list);

    double x, y;
    for (int i = 0; i < _n_cam_all; i ++)
    {
        x = obj3d._tr2d_list[i]._pt_center[0];
        y = obj3d._tr2d_list[i]._pt_center[1];

        if (x < 0 || x >= _cam_list.cam_list[i].getNCol() ||
            y < 0 || y >= _cam_list.cam_list[i].getNRow())
        {
            n_outrange ++;

            if (_n_cam_all - n_outrange < 2)
            {
                return false;
            }
        }
    }

    return true;
}


template<class T3D>
int STB<T3D>::linkShortTrack (Track<T3D> const& track, std::vector<T3D> const& obj3d_list, int n_iter)
{
    int n_la = _long_track_active.size();
    Pt3D est, vel;
    double w, beta2;
    int obj3d_id = UNLINKED;

    double tot_weight = 0;
    std::vector<double> weight;
    int n_nearTracks = 0;
    std::vector<Pt3D> disp;

    std::vector<int> is_around(n_la, 0);
    std::vector<double> dsqr(n_la, 1e8);
    double rsqr, shift;
    
    for (int step = 0; step < n_iter; step ++)
    {
        rsqr = std::pow(std::pow(1.1, step) * 3 * _r_trackSearch, 2);
        shift = std::pow(1.1, step) * _r_objSearch;
        
        vel *= 0;
        beta2 = 1 / std::pow(std::pow(1.1, step) * 0.1 * _r_trackSearch, 2); // typical length scale
        
        // calculate the predictive vel. field as an avg. of particle vel from neighbouring tracks
        // identify the neighbour tr (using 3*avg interpt dist.) and get their pt vel
        int len = 0;
        for (int j = 0; j < n_la; j ++)
        {
            if (is_around[j])
            {
                continue;
            }

            len = _long_track_active[j]._obj3d_list.size();

            if (step == 0)
            {
                dsqr[j] = myMATH::dist2(
                    _long_track_active[j]._obj3d_list[len-2]._pt_center, 
                    track._obj3d_list.back()._pt_center
                ); 
            }

            if (dsqr[j] < rsqr)
            {
                is_around[j] = 1;
                
                w = 1/(1 + dsqr[j]*beta2);
                // w = std::sqrt(dsqr[j]);
                tot_weight += w;
                weight.push_back(w);
                
                disp.push_back(
                    _long_track_active[j]._obj3d_list[len-1]._pt_center - _long_track_active[j]._obj3d_list[len-2]._pt_center
                );
            }
        }

        n_nearTracks = weight.size();
        if (n_nearTracks > 0)
        {
            for (int j = 0; j < n_nearTracks; j ++)
            {
                w = weight[j] / tot_weight;              
                vel[0] += w * disp[j][0];
                vel[1] += w * disp[j][1];
                vel[2] += w * disp[j][2];
            }

            est = track._obj3d_list.back()._pt_center + vel;
            obj3d_id = findNN(obj3d_list, est, _r_objSearch);
        }
        else
        {
            // if no neighbouring tracks are identified
            est = track._obj3d_list.back()._pt_center;
            obj3d_id = findNN(obj3d_list, est, shift);
        }

        if (obj3d_id != UNLINKED)
        {
            break;
        }
    }

    return obj3d_id;
}


template<class T3D>
bool STB<T3D>::checkLinearFit(Track<T3D> const& track)
{
    int len = track._obj3d_list.size();
    int n_pts = _n_initPhase > 4 ? 4 : _n_initPhase;

    if (len < n_pts)
    {
        return false;
    }

    Matrix<double> coeff(3, 2, 0);
    Matrix<double> x_mat(n_pts, 2, 0);
    Matrix<double> kernel(2, n_pts, 0);
    Matrix<double> y_mat(n_pts, 1, 0);
    Matrix<double> temp(2, 1, 0);

    // prepare x_mat
    for (int i = 0; i < n_pts; i ++)
    {
        x_mat(i, 0) = 1;
        x_mat(i, 1) = i;
    }

    // calculate the kernel matrix
    kernel = myMATH::inverse(x_mat.transpose() * x_mat) * x_mat.transpose();

    // calculate coeff at each direction and the estimate point
    Pt3D est;
    for (int i = 0; i < 3; i ++)
    {
        for (int j = 0; j < n_pts; j ++)
        {
            y_mat(j, 0) = track._obj3d_list[len-n_pts+j]._pt_center[i];
        }

        temp = kernel * y_mat;
        coeff(i, 0) = temp(0, 0);
        coeff(i, 1) = temp(1, 0);
        est[i] = coeff(i, 0) + (n_pts-1) * coeff(i, 1);
    }

    // calculate the residue
    double res = myMATH::dist(track._obj3d_list[len-1]._pt_center, est);
    double err_max = MAX_ERR_LINEARFIT;
    double err_min = 0.6 * err_max;

    // std::cout << "residue: " << res2 << std::endl;

    if (res > err_max)
    {
        return false;
    }
    else if (res > err_min)
    {
        // For those tracks that are at the edge of being deleted (res>0.6*err), check the whole track to see whether its error is large also
        double res_tot = res;
        for (int i = 0; i < n_pts-1; i ++)
        {
            est[0] = coeff(0,0) + coeff(0,1) * i;
            est[1] = coeff(1,0) + coeff(1,1) * i;
            est[2] = coeff(2,0) + coeff(2,1) * i;

            res = myMATH::dist(track._obj3d_list[len-n_pts+i]._pt_center, est); 

            if (res > err_min)
            {
                return false;
            }

            res_tot += res;
        }
        if (res_tot > n_pts * err_min)
        {
            return false;
        }
    }

    return true;
}


template<class T3D>
void STB<T3D>::loadTracks (std::string const& file, TrackStatusID status)
{
    std::ifstream input(file, std::ios::in);

    if (!input.is_open())
    {
        std::cerr << "STB<T3D>::loadTracks error at line" << __LINE__ << ":\n"
                  << "Cannot open file " << file << std::endl;
        throw error_io;
    }

    int track_id_prev = -1;
    std::vector<double> data(5,0); // track_id, frame_id, x, y, z  
    T3D obj3d;
    if (typeid(T3D) == typeid(Tracer3D))
    {
        obj3d._r2d_px = _obj_param[2];
    }

    std::deque<Track<T3D>> track_list;
    std::string line;
    std::getline(input, line); // skip the first line
    while (std::getline(input, line))
    {
        // Trim leading and trailing whitespace (optional)
        line.erase(0, line.find_first_not_of(" \t\r\n")); // Trim leading whitespace
        line.erase(line.find_last_not_of(" \t\r\n") + 1); // Trim trailing whitespace
        if (line.empty())
        {
            continue;
        }
        
        std::istringstream ss(line);

        for (int i = 0; i < 5; i ++)
        {
            ss >> data[i];
            ss.ignore();
        }
        
        obj3d._pt_center[0] = data[2];
        obj3d._pt_center[1] = data[3];
        obj3d._pt_center[2] = data[4];

        if (data[0] != track_id_prev)
        {
            Track<T3D> track(obj3d, data[1], _track_pred_param);
            track_list.push_back(track);
            track_id_prev = data[0];
        }
        else
        {
            track_list.back().addNext(obj3d, data[1]);
        }
    }

    input.close();


    switch (status)
    {
        case (LONG_ACTIVE):
            _long_track_active = track_list;
            break;
        case (LONG_INACTIVE):
            _long_track_inactive = track_list;
            break;
        case (SHORT_ACTIVE):
            _short_track_active = track_list;
            break;
        case (EXIT):
            _exit_track = track_list;
            break;
        default:
            std::cerr << "STB<T3D>::loadTracks error at line" 
                      << __LINE__ << ":\n"
                      << "Invalid track status: " 
                      << status << std::endl;
            throw error_type;
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

    output << "TrackID,FrameID,WorldX,WorldY,WorldZ,Error,Ncam";
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