#ifndef IPR_CPP
#define IPR_CPP

#include "IPR.h"


void IPR::runIPR(
    std::vector<Tracer3D>& tr3d_list_all, 
    std::vector<double> const& tr2d_properties, 
    OTF const& otf, 
    bool is_reduced, int n_reduced) 
{
    clock_t t_tot_start, t_tot_end;
    t_tot_start = clock();

    // reset cam id list
    resetCamIDList ();

    // Clear output
    tr3d_list_all.clear();

    // Initialize stereo match parameters
    StereoMatchParam match_param;
    match_param.tor_2d = _param.tol_2d;
    match_param.tor_3d = _param.tol_3d;
    match_param.n_thread = _param.n_thread;
    match_param.check_id = _param.check_id < _n_cam_all ? _param.check_id : _n_cam_all;
    match_param.check_radius = _param.check_radius;
    match_param.is_delete_ghost = true;
    match_param.is_update_inner_var = false;
    StereoMatch stereo_match(match_param, _cam_list);

    // Initialize shake parameters
    Shake s (
        _cam_list, 
        _param.shake_width, 
        _param.ghost_threshold, 
        _param.n_loop_shake, 
        _param.n_thread
    ); // gradient descent


    // Start IPR loop
    ObjectFinder2D objfinder;
    for (int loop = 0; loop < _param.n_loop_ipr; loop ++)
    {
        // Step 1: identify 2D position from image
        std::vector<std::vector<Tracer2D>> tr2d_list_all;
        std::cout << "\tNumber of found tracers in each camera: ";
        for (int i = 0; i < _n_cam_all; i ++)
        {
            std::vector<Tracer2D> tr2d_list;
            objfinder.findObject2D(tr2d_list, _imgRes_list[i], tr2d_properties);
            tr2d_list_all.push_back(tr2d_list);

            std::cout << tr2d_list.size() << ",";
        }
        std::cout << std::endl;


        // Step 2: stereomatching
        std::vector<Tracer3D> tr3d_list;
        clock_t t_start, t_end;
        t_start = clock();
        stereo_match.match(tr3d_list, tr2d_list_all);
        t_end = clock();
        std::cout << "\tMatching time = " 
                  << (double) (t_end - t_start)/CLOCKS_PER_SEC
                  << " [s]" << std::endl;


        // Step 3: Shake
        if (tr3d_list.size() == 0)
        {
            break;
        }
        t_start = clock();
        s.runShake(tr3d_list, otf, _imgRes_list, _param.tri_only); // 0.25 vox, 1 vox = 0.04 mm
        t_end = clock();
        std::cout << "\tShake time = " 
                  << double(t_end-t_start)/CLOCKS_PER_SEC 
                  << " [s]" << std::endl;


        // Save tracer info after shaking
        tr3d_list_all.insert(tr3d_list_all.end(), tr3d_list.begin(), tr3d_list.end());

        // Update imgRes_list
        _imgRes_list = s._imgRes_list;

        std::cout << "\tIPR step " << loop << ": find " << tr3d_list_all.size() << " particles. " << std::endl;
    }


    // Reduced cam
    if (is_reduced)
    {
        int n_rest = _n_cam_all - n_reduced;
        if (_n_cam_all - n_reduced < 2 || n_reduced < 0)
        {
            std::cerr << "IPR::runIPR error at line " << __LINE__ << ": \n"
                      << "_n_cam_all = " << _n_cam_all 
                      << ", n_reduced = " << n_reduced << "\n"
                      << "_n_cam_all - n_reduced = " << n_rest << "\n"
                      << "It is out of range."
                      << std::endl;
            throw error_range;
        }

        std::cout << "\tStart reduced camera loop!" << std::endl;

        // Generate all possible cam_id
        std::deque<std::vector<int>> cam_id_all;
        std::vector<int> cam_id;
        createCamID (cam_id_all, cam_id, 0, n_rest);

        for (int i = 0; i < cam_id_all.size(); i ++)
        {
            reducedCamLoop(tr3d_list_all, tr2d_properties, otf, cam_id_all[i], _param.n_loop_ipr_reduced);
        }
    }

    // reset cam id list
    resetCamIDList ();

    t_tot_end = clock();
    std::cout << "IPR Finish! " << "Find " << tr3d_list_all.size() << " particles. " << " Total time = " << (t_tot_end-t_tot_start)/CLOCKS_PER_SEC << " [s]" << std::endl;
}


void IPR::reducedCamLoop(std::vector<Tracer3D>& tr3d_list_all, std::vector<double> const& tr2d_properties, 
    OTF const& otf, std::vector<int> const& cam_id, int n_loop)
{
    _cam_list.useid_list = cam_id;

    // Clear output
    tr3d_list_all.clear();

    // Initialize stereo match parameters
    StereoMatchParam match_param;
    match_param.tor_2d = _param.tol_2d;
    match_param.tor_3d = _param.tol_3d;
    match_param.n_thread = _param.n_thread;
    match_param.check_id = _param.check_id < _n_cam_all ? _param.check_id : _n_cam_all;
    match_param.check_radius = _param.check_radius;
    match_param.is_delete_ghost = true;
    match_param.is_update_inner_var = false;
    StereoMatch stereo_match(match_param, _cam_list);

    // Initialize shake parameters
    Shake s (
        _cam_list, 
        _param.shake_width, 
        _param.ghost_threshold, 
        _param.n_loop_shake, 
        _param.n_thread
    ); // gradient descent

    // Start IPR loop
    ObjectFinder2D objfinder;
    for (int loop = 0; loop < n_loop; loop ++)
    {
        // Step 1: identify 2D position from image
        std::vector<std::vector<Tracer2D>> tr2d_list_all;
        std::cout << "\tNumber of found tracers in each camera: ";
        for (int i = 0; i < _n_cam_all; i ++)
        {
            std::vector<Tracer2D> tr2d_list;
            objfinder.findObject2D(tr2d_list, _imgRes_list[i], tr2d_properties);
            tr2d_list_all.push_back(tr2d_list);

            std::cout << tr2d_list.size() << ",";
        }
        std::cout << std::endl;


        // Step 2: stereomatching
        std::vector<Tracer3D> tr3d_list;
        clock_t t_start, t_end;
        t_start = clock();
        stereo_match.match(tr3d_list, tr2d_list_all);
        t_end = clock();
        std::cout << "\tMatching time = " 
                  << (double) (t_end - t_start)/CLOCKS_PER_SEC
                  << " [s]" << std::endl;


        // Step 3: Shake
        if (tr3d_list.size() == 0)
        {
            break;
        }
        t_start = clock();
        s.runShake(tr3d_list, otf, _imgRes_list, _param.tri_only); // 0.25 vox, 1 vox = 0.04 mm
        t_end = clock();
        std::cout << "\tShake time = " 
                  << double(t_end-t_start)/CLOCKS_PER_SEC 
                  << " [s]" << std::endl;


        // Save tracer info after shaking
        tr3d_list_all.insert(tr3d_list_all.end(), tr3d_list.begin(), tr3d_list.end());

        // Update imgRes_list
        _imgRes_list = s._imgRes_list;

        std::cout << "\tIPR step " << loop << ": find " << tr3d_list_all.size() << " particles. " << std::endl;
    }
    
}


void IPR::createCamID (std::deque<std::vector<int>>& cam_id_all, std::vector<int> cam_id, int id, int n_rest)
{
    int n_cam_id = cam_id.size();
    if (n_cam_id < n_rest && id < _n_cam_all-(n_rest-n_cam_id)+1)
    {
        for (int i = id; i < _n_cam_all-(n_rest-n_cam_id)+1; i ++)
        {
            cam_id.push_back(i);
            createCamID(cam_id_all, cam_id, i+1, n_rest);
            cam_id.pop_back();
        }
    }
    else
    {
        if (n_cam_id == n_rest)
        {
            cam_id_all.push_back(cam_id);
        }
        else
        {
            std::cerr << "IPR::createCamID problem!" 
                      << "cam_id.size=" << cam_id.size()
                      << ", n_rest=" << n_rest 
                      << std::endl;
            throw error_range;
        }
    }
}


void IPR::resetCamIDList ()
{
    _cam_list.useid_list.resize(_n_cam_all);
    for (int i = 0; i < _n_cam_all; i ++)
    {
        _cam_list.useid_list[i] = i;
    }
}


#endif