#ifndef IPR_CPP
#define IPR_CPP

#include "IPR.h"


void IPR::runIPR(
    std::vector<Tracer3D>& tr3d_list_all, 
    std::vector<double> const& tr2d_properties, 
    OTF const& otf, 
    int n_reduced) 
{
    clock_t t_tot_start, t_tot_end;
    t_tot_start = clock();

    // clear output
    tr3d_list_all.clear();

    // reset cam id list
    resetCamIDList ();

    // Initialize stereo match parameters
    SMParam match_param;
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
        _param.tol_3d * 3,
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
            std::cout << tr2d_list.size();

            // if tr2d_list is too large, randomly select some tracers
            int seed = 1234;
            if (tr2d_list.size() > _param.n_obj2d_max)
            {
                std::shuffle(tr2d_list.begin(), tr2d_list.end(), std::default_random_engine(seed));
                tr2d_list.erase(tr2d_list.begin()+_param.n_obj2d_max, tr2d_list.end());
                std::cout << "(" << _param.n_obj2d_max << ")";
            }

            tr2d_list_all.push_back(tr2d_list);

            std::cout << ",";
            
            if (tr2d_list.size() == 0)
            {
                std::cout << "\n\tQuit IPR: No tracer found in camera " << i << std::endl;
                return;
            }
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
        for (int i = tr3d_list.size()-1; i >= 0; i --)
        {
            if (s._is_ghost[i])
            {
                tr3d_list.erase(tr3d_list.begin()+i);
            }
        }
        tr3d_list_all.insert(tr3d_list_all.end(), tr3d_list.begin(), tr3d_list.end());

        // Update imgRes_list
        // Note: s._imgRes_list in shake is the same size as n_cam_use, but _imgRes_list in IPR is the same size as _n_cam_all.
        // Since there is no reduced camera here, s._imgRes_list is the same as _imgRes_list.
        _imgRes_list = s._imgRes_list;

        std::cout << "  IPR step " << loop << ": find " << tr3d_list_all.size() << " particles. " << std::endl;

        if (tr3d_list.size() == 0)
        {
            break;
        }
    }


    // Reduced cam
    if (n_reduced > 0)
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

        std::cout << "Start reduced camera loop!" << std::endl;

        for (int n_reduced_cur = 1; n_reduced_cur <= n_reduced; n_reduced_cur ++)
        {
            int n_rest_cur = _n_cam_all - n_reduced_cur;
            std::cout << "Reduce: " << n_reduced_cur << " cam:" << std::endl;

            // Generate all possible cam_id
            std::deque<std::vector<int>> cam_id_all;
            std::vector<int> cam_id;
            createCamID (cam_id_all, cam_id, 0, n_rest_cur);
            std::reverse(cam_id_all.begin(), cam_id_all.end());

            reducedCamLoop(tr3d_list_all, tr2d_properties, otf, cam_id_all, _param.n_loop_ipr_reduced);
            
        } 
    }

    // reset cam id list
    resetCamIDList ();

    t_tot_end = clock();
    std::cout << "IPR Finish! " << "Find " << tr3d_list_all.size() << " particles." << " Total time = " << (double)(t_tot_end-t_tot_start)/CLOCKS_PER_SEC << " [s]" << std::endl;
}


void IPR::reducedCamLoop(std::vector<Tracer3D>& tr3d_list_all, std::vector<double> const& tr2d_properties, 
    OTF const& otf, std::deque<std::vector<int>> const& cam_id_all, int n_loop)
{
    // Start IPR loop
    ObjectFinder2D objfinder;
    int cam_id;
    for (int loop = 0; loop < n_loop; loop ++)
    {
        for (int combine_id = 0; combine_id < cam_id_all.size(); combine_id ++)
        {   
            _cam_list.useid_list = cam_id_all[combine_id];
            int n_cam_use = _cam_list.useid_list.size();

            // Initialize stereo match parameters
            SMParam match_param;
            match_param.tor_2d = _param.tol_2d;
            match_param.tor_3d = _param.tol_3d;
            match_param.n_thread = _param.n_thread;
            match_param.check_id = _param.check_id < n_cam_use ? _param.check_id : n_cam_use;
            match_param.check_radius = _param.check_radius;
            match_param.is_delete_ghost = true;
            match_param.is_update_inner_var = false;
            StereoMatch stereo_match(match_param, _cam_list);

            // Initialize shake parameters
            Shake s (
                _cam_list, 
                _param.shake_width, 
                _param.tol_3d * 3,
                _param.ghost_threshold, 
                _param.n_loop_shake, 
                _param.n_thread
            ); // gradient descent

            // Step 1: identify 2D position from image
            std::vector<std::vector<Tracer2D>> tr2d_list_all;
            std::cout << "\tNumber of found tracers in each camera: ";
            for (int i = 0; i < n_cam_use; i ++)
            {
                cam_id = _cam_list.useid_list[i];

                std::vector<Tracer2D> tr2d_list;
                objfinder.findObject2D(tr2d_list, _imgRes_list[cam_id], tr2d_properties);

                std::cout << tr2d_list.size();

                // if tr2d_list is too large, randomly select some tracers
                int seed = 123;
                if (tr2d_list.size() > _param.n_obj2d_max)
                {
                    std::shuffle(tr2d_list.begin(), tr2d_list.end(), std::default_random_engine(seed));
                    tr2d_list.erase(tr2d_list.begin()+_param.n_obj2d_max, tr2d_list.end());

                    std::cout << "(" << _param.n_obj2d_max << ")";
                }

                tr2d_list_all.push_back(tr2d_list);

                std::cout << ",";
                
                if (tr2d_list.size() == 0)
                {
                    std::cout << "\n\tQuit IPR Reduced camera: No tracer found in camera " << cam_id << std::endl;
                    return;
                }
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
            for (int i = tr3d_list.size()-1; i >= 0; i --)
            {
                if (s._is_ghost[i])
                {
                    tr3d_list.erase(tr3d_list.begin()+i);
                }
            }
            tr3d_list_all.insert(tr3d_list_all.end(), tr3d_list.begin(), tr3d_list.end());


            // Update imgRes_list
            // Note: s._imgRes_list in shake is the same size as n_cam_use, but _imgRes_list in IPR is the same size as _n_cam_all
            // Since there is reduced camera here, s._imgRes_list is not the same as _imgRes_list
            for (int i = 0; i < n_cam_use; i ++)
            {
                cam_id = _cam_list.useid_list[i];
                _imgRes_list[cam_id] = s._imgRes_list[i];
            }


            std::cout << "\tIPR reduced camera step " << loop << ": find " << tr3d_list_all.size() << " particles. " << std::endl;

        }
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