#ifndef STEREOMATCH_HPP
#define STEREOMATCH_HPP

#include "StereoMatch.h"

StereoMatch::StereoMatch(SMParam const& param, CamList const& cam_list) 
    : _param(param), _cam_list(cam_list), _n_cam_use(cam_list.useid_list.size()) 
{
    if (_n_cam_use < 2)
    {
        std::cerr << "StereoMatch::StereoMatch error at line " << __LINE__ << ": \n"
                  << "Need at least 2 cameras for matching" 
                  << std::endl;
        throw error_size;
    }
    else if (param.check_id > _n_cam_use)
    {
        std::cerr << "StereoMatch::StereoMatch error at line " << __LINE__ << ": \n"
                  << "param.check_id = " << param.check_id 
                  << " > _n_cam_use = " << _n_cam_use 
                  << std::endl;
        throw error_size;
    }
}

void StereoMatch::clearAll()
{
    _objID_map_list.clear();
    _error_list.clear();
    _objID_match_list.clear();
    _n_before_del = 0;
}

template<class T3D, class T2D>
void StereoMatch::match(std::vector<T3D>& obj3d_list, std::vector<std::vector<T2D>> const& obj2d_list)
{
    // clear all the lists
    clearAll();
    obj3d_list.clear();

    if (typeid(T3D)==typeid(Tracer3D) && typeid(T2D)==typeid(Tracer2D))
    {
        tracerMatch(obj2d_list);

        if (_param.is_delete_ghost)
        {   
            // removeGhostTracer(obj3d_list, obj2d_list);
            removeGhostTracerTest(obj3d_list, obj2d_list);
        }
        else
        {
            fillTracerInfo(obj3d_list, obj2d_list);
        }
    }
    else
    {
        std::cerr << "StereoMatch::match error at line " << __LINE__ << ": \n"
                  << "The type of object is not supported." 
                  << std::endl;
        throw error_type;
    }
}

template<class T3D>
void StereoMatch::saveObjInfo (std::string path, std::vector<T3D> const& obj3d_list)
{
    if (typeid(T3D)==typeid(Tracer3D))
    {
        saveTracerInfo(path, obj3d_list);
    }
    else
    {
        std::cerr << "StereoMatch::saveObjInfo error at line " << __LINE__ << ": \n"
                  << "The type of object is not supported." 
                  << std::endl;
        throw error_type;
    }
}

// save obj ID match list
void StereoMatch::saveObjIDMatchList (std::string path)
{
    std::ofstream file;
    file.open(path, std::ios::out);

    for (int i = 0; i < _objID_match_list.size(); i ++)
    {
        for (int j = 0; j < _n_cam_use-1; j ++)
        {
            file << _objID_match_list[i][j] << ",";
        }
        file << _objID_match_list[i][_n_cam_use-1] << "\n";
    }

    file.close();
}

// create object ID map
template<class T2D>
void StereoMatch::createObjIDMap (std::vector<std::vector<T2D>> const& obj2d_list)
{
    int row_id, col_id;
    int cam_id;
    for (int i = 0; i < _n_cam_use; i ++)
    {
        cam_id = _cam_list.useid_list[i];

        ObjIDMap objID_map;
        objID_map.config(
            _cam_list.cam_list[cam_id].getNRow(), 
            _cam_list.cam_list[cam_id].getNCol()
        );
        
        for (int j = 0; j < obj2d_list[i].size(); j ++)
        {
            row_id = obj2d_list[i][j]._pt_center[1]; // img_y
            col_id = obj2d_list[i][j]._pt_center[0]; // img_x

            if (objID_map(row_id, col_id)[0] == -1)
            {
                objID_map(row_id, col_id)[0] = j;
            }
            else
            {
                objID_map(row_id, col_id).push_back(j);
            }
        }

        _objID_map_list.push_back(objID_map);
    }
}


//              //
// Tracer match //
//              //

// get list of matched tracer_id list
void StereoMatch::tracerMatch (std::vector<std::vector<Tracer2D>> const& tr2d_list)
{
    if (_n_cam_use < 2)
    {
        std::cerr << "StereoMatch::tracerMatch error at line " << __LINE__ << ": \n"
                  << "Need at least 2 cameras for matching" 
                  << std::endl;
        throw error_size;
    }
    if (int(tr2d_list.size()) != _n_cam_use)
    {
        std::cerr << "StereoMatch::tracerMatch error at line " << __LINE__ << ":\n"
                  << "The size of tracer list is: " 
                  << tr2d_list.size() 
                  << ", is different from the size of camera list: "
                  << _n_cam_use << "."
                  << std::endl;
        throw error_size;
    }

    createObjIDMap (tr2d_list);

    if (_param.n_thread != 0)
    {
        omp_set_num_threads(_param.n_thread);
    }
    #pragma omp parallel
    {
        // for 1st used camera, draw a line of sight through each particle on its image plane.
        // project these lines of sight onto the image planes of 2nd camera.
        // particles within mindist_2D of these lines are candidate matches between first 2 cams.
        // then project 2 line of sights from each particle pair of 1st & 2nd cam onto 3rd cam.
        // particles within torlerance are candidate matches from 3rd cam.
        // repeat similarly for subsequent cams
        #pragma omp for 
        for (int tr_id = 0; tr_id < tr2d_list[0].size(); tr_id ++)
        {
            std::deque<std::vector<int>> trID_match_list; // match list for the particle i in the first camera

            std::vector<int> trID_match; 
            trID_match.push_back(tr_id);

            std::deque<double> error_list;
            
            findTracerMatch (
                1,                
                trID_match,
                trID_match_list,
                error_list,
                tr2d_list
            );

            #pragma omp critical
            {   
                if (error_list.size() != trID_match_list.size())
                {
                    std::cerr << "StereoMatch::tracerMatch error at line " << __LINE__ << ": \n" 
                              << "error_list.size = " << error_list.size()
                              << "tracer_id_match_list.size = " << trID_match_list.size() << "\n";
                    throw error_size;
                }
                for (int i = 0; i < trID_match_list.size(); i ++) 
                {      
                    if (trID_match_list[i].size() != _n_cam_use)
                    {
                        continue;
                    }

                    int n_error_size_pre = _error_list.size();
                    int n_match_size_pre = _objID_match_list.size();

                    _error_list.push_back(error_list[i]);
                    _objID_match_list.push_back(trID_match_list[i]);  
                    
                    int n_error_size = _error_list.size();
                    int n_match_size = _objID_match_list.size();
                    if (n_error_size == n_error_size_pre || n_match_size == n_match_size_pre)
                    {
                        std::cerr << "StereoMatch::tracerMatch error at line " << __LINE__ << ":\n"
                                  << "_error_list.size != _object_id_match_list.size; "
                                  << "_error_list.size = " << n_error_size
                                  << "," << n_error_size_pre << ", _object_id_match_list.size = " << n_match_size << "," << n_match_size_pre
                                  << "\n";
                        throw error_size;
                    }
                }
            }
        }
    }
    
    if (_error_list.size() != _objID_match_list.size())
    {
        std::cerr << "StereoMatch::tracerMatch error at line " << __LINE__ << ":\n"
                  << "_error_list.size != _objID_match_list.size; "
                  << "_error_list.size = " << _error_list.size()
                  << ", _object_id_match_list.size = " << _objID_match_list.size()
                  << "\n";
        throw error_size;
    }

    _n_before_del = _objID_match_list.size();
    std::cout << "\tFinish stereomatching: " 
              << "n_before_del = " << _n_before_del << "."
              << std::endl;
}

// remove ghost tracer
void StereoMatch::removeGhostTracer (std::vector<Tracer3D>& tr3d_list, std::vector<std::vector<Tracer2D>> const& tr2d_list)
{
    int n_match = _objID_match_list.size();
    if (n_match == 0)
    {
        return;
    }

    std::vector<int> match_score(n_match, 0);

    // update match_score  
    int n_tr2d, tr2d_id, opt_matchID;
    for (int i = 0; i < _n_cam_use; i ++)
    {
        n_tr2d = tr2d_list[i].size();

        // optMatchID_map: 
        //  save the corresponding match id with smallest error for each tracer2d 
        std::vector<int> optMatchID_map(n_tr2d, -1);
        
        // update optMatchID_map
        for (int j = 0; j < n_match; j ++)
        {
            tr2d_id = _objID_match_list[j][i];
            opt_matchID = optMatchID_map[tr2d_id];

            if (opt_matchID == -1) 
            {
                optMatchID_map[tr2d_id] = j;
            }
            else 
            {
                if (_error_list[opt_matchID] > _error_list[j])
                {
                    // replace by a match id with smaller error 
                    optMatchID_map[tr2d_id] = j;
                }
            }
        }


        // update the score
        for (int j = 0; j < n_tr2d; j ++)
        {
            if (optMatchID_map[j] != -1)
            {
                opt_matchID = optMatchID_map[j];
                match_score[opt_matchID] += 1;
            }
        }

    }

    // sort match_score from large to small
    //  small to large 
    std::vector<int> matchID_list(n_match);
    myMATH::sortID(matchID_list, match_score);

    // map tr2d_id to match_id with smallest error
    std::vector<std::vector<int>> tr2dID_matchID_map;
    for (int i = 0; i < _n_cam_use; i ++)
    {
        n_tr2d = tr2d_list[i].size();
        tr2dID_matchID_map.push_back(std::vector<int>(n_tr2d, -1));
    }

    // select the optimal match
    //  from the largest score to the smallest score
    std::vector<bool> is_select(n_match, false);
    int match_id, map_matchID, n_equal;
    double sum_error;
    bool is_rank_small;
    for (int i = n_match-1; i > -1; i --)
    {
        match_id = matchID_list[i];  
        n_equal = 0;
        sum_error = 0.0;
        is_rank_small = false;

        for (int j = 0; j < _n_cam_use; j ++)
        {
            tr2d_id = _objID_match_list[match_id][j];
            map_matchID = tr2dID_matchID_map[j][tr2d_id];

            if (map_matchID == -1)
            {
                // when there is no match for the checked particle
                continue;          
            }
            else if ((match_score[match_id] < match_score[map_matchID]) && (!is_rank_small)) 
            {
                is_rank_small = true;
                break;
            }
            else if (match_score[match_id] == match_score[map_matchID])
            {
                // when the preference # is the same as the current preference # of those filled matches,
                //  then we should calculate the average of the preference error of those filled matches
                //  and compare it with the current preference error
                n_equal ++;
                sum_error += _error_list[map_matchID];
            }
        }

        if ( !is_rank_small &&
            (( n_equal > 0 && _error_list[match_id] < sum_error/n_equal ) || 
              n_equal == 0 ))
        {
            // 1) if there is same preference number (for certain cam),
            // and the error is less than the averaged of those filled match, 
            // then replace all the old match with the new one
            // 2) no repeating

            for (int j = 0; j < _n_cam_use; j ++)
            {
                tr2d_id = _objID_match_list[match_id][j];
                map_matchID = tr2dID_matchID_map[j][tr2d_id];

                if (map_matchID != -1)
                {
                    is_select[map_matchID] = false;
                    for (int k = 0; k < _n_cam_use; k ++)
                    {
                        int map_matchID_orig = _objID_match_list[map_matchID][k];
                        tr2dID_matchID_map[k][map_matchID_orig] = -1;
                    }
                }

                tr2dID_matchID_map[j][tr2d_id] = match_id;
            }

            is_select[match_id] = true;
        }
        
    }
    
    // get new match list 
    Tracer3D tr3d;
    tr3d._camid_list = _cam_list.useid_list;
    tr3d._n_2d = _n_cam_use;
    tr3d._tr2d_list.resize(_n_cam_use);
    std::vector<Line3D> sight3D_list(_n_cam_use);

    int cam_id;
    std::vector<std::vector<int>> objID_match_list_new;
    std::vector<double> error_list_new;
    for (int i = 0; i < n_match; i ++)
    {
        if (is_select[i])
        {
            for (int id = 0; id < _n_cam_use; id ++)
            {
                tr2d_id = _objID_match_list[i][id];
                cam_id = _cam_list.useid_list[id];

                tr3d._tr2d_list[id]._pt_center = tr2d_list[id][tr2d_id]._pt_center;
                tr3d._tr2d_list[id]._r_px = tr2d_list[id][tr2d_id]._r_px;
                sight3D_list[id] = _cam_list.cam_list[cam_id].lineOfSight(tr3d._tr2d_list[id]._pt_center);
            }

            tr3d._r2d_px = tr3d._tr2d_list[0]._r_px;
            myMATH::triangulation(tr3d._pt_center, tr3d._error, sight3D_list);

            tr3d_list.push_back(tr3d);

            if (_param.is_update_inner_var)
            {
                objID_match_list_new.push_back(_objID_match_list[i]);
                error_list_new.push_back(_error_list[i]);
            }
        }
    }
    if (_param.is_update_inner_var)
    {
        _objID_match_list = objID_match_list_new;
        _error_list = error_list_new;
    }

    // print info
    _n_after_del = tr3d_list.size();
    _n_del = _n_before_del - _n_after_del;
    std::cout << "\tFinish deleting gohst match: "
              << "n_del = " << _n_del << ", "
              << "n_after_del = " << _n_after_del << "."
              << std::endl;
}

// remove ghost tracer
void StereoMatch::removeGhostTracerTest (std::vector<Tracer3D>& tr3d_list, std::vector<std::vector<Tracer2D>> const& tr2d_list)
{
    int n_match = _objID_match_list.size();
    if (n_match == 0)
    {
        return;
    }

    std::vector<int> match_score(n_match, 0);

    // update match_score  
    int n_tr2d, tr2d_id, opt_matchID;
    for (int i = 0; i < _n_cam_use; i ++)
    {
        n_tr2d = tr2d_list[i].size();

        // optMatchID_map: 
        //  save the corresponding match id with smallest error for each tracer2d 
        std::vector<int> optMatchID_map(n_tr2d, -1);
        
        // update optMatchID_map
        for (int j = 0; j < n_match; j ++)
        {
            tr2d_id = _objID_match_list[j][i];
            opt_matchID = optMatchID_map[tr2d_id];

            if (opt_matchID == -1) 
            {
                optMatchID_map[tr2d_id] = j;
            }
            else 
            {
                if (_error_list[opt_matchID] > _error_list[j])
                {
                    // replace by a match id with smaller error 
                    optMatchID_map[tr2d_id] = j;
                }
            }
        }


        // update the score
        for (int j = 0; j < n_tr2d; j ++)
        {
            if (optMatchID_map[j] != -1)
            {
                opt_matchID = optMatchID_map[j];
                match_score[opt_matchID] += 1;
            }
        }

    }

    // sort match_score and error_list
    //  match_score: small to large 
    //  error_list: large to small
    std::vector<int> matchID_list(n_match);
    for (int i = 0; i < n_match; i ++)
    {
        matchID_list[i] = i;
    }
    std::vector<double>& error_list(_error_list);
    auto comparator = [match_score, error_list](size_t i, size_t j) {
        if (match_score[i] < match_score[j])
        {
            return true;
        }
        else if (match_score[i] == match_score[j])
        {
            return error_list[i] > error_list[j];
        }
        else
        {
            return false;
        }
    };
    std::sort(matchID_list.begin(), matchID_list.end(), comparator);

    // create map for record whether a tr2d is used
    std::vector<std::vector<bool>> is_tr2d_use;
    for (int i = 0; i < _n_cam_use; i ++)
    {
        n_tr2d = tr2d_list[i].size();
        is_tr2d_use.push_back(std::vector<bool>(n_tr2d, false));
    }

    // select the optimal match based on the sortID
    //  from the highest score and smallest error to the lowest score and largest error
    //  from the last one to the first one
    std::vector<bool> is_select(n_match, true);
    int match_id;
    bool is_use;
    for (int i = n_match-1; i > -1; i --)
    {
        match_id = matchID_list[i];

        for (int j = 0; j < _n_cam_use; j ++)
        {
            tr2d_id = _objID_match_list[match_id][j];
            is_use = is_tr2d_use[j][tr2d_id];
            if (is_use)
            {
                is_select[match_id] = false;
                break;
            }
        }

        if (is_select[match_id])
        {
            for (int j = 0; j < _n_cam_use; j ++)
            {
                tr2d_id = _objID_match_list[match_id][j];
                is_tr2d_use[j][tr2d_id] = true;
            }
        }
    }
    
    // get new match list 
    Tracer3D tr3d;
    tr3d._camid_list = _cam_list.useid_list;
    tr3d._n_2d = _n_cam_use;
    tr3d._tr2d_list.resize(_n_cam_use);
    std::vector<Line3D> sight3D_list(_n_cam_use);

    int cam_id;
    std::vector<std::vector<int>> objID_match_list_new;
    std::vector<double> error_list_new;
    for (int i = 0; i < n_match; i ++)
    {
        if (is_select[i])
        {
            for (int id = 0; id < _n_cam_use; id ++)
            {
                tr2d_id = _objID_match_list[i][id];
                cam_id = _cam_list.useid_list[id];

                tr3d._tr2d_list[id]._pt_center = tr2d_list[id][tr2d_id]._pt_center;
                tr3d._tr2d_list[id]._r_px = tr2d_list[id][tr2d_id]._r_px;
                sight3D_list[id] = _cam_list.cam_list[cam_id].lineOfSight(tr3d._tr2d_list[id]._pt_center);
            }

            tr3d._r2d_px = tr3d._tr2d_list[0]._r_px;
            myMATH::triangulation(tr3d._pt_center, tr3d._error, sight3D_list);

            tr3d_list.push_back(tr3d);

            if (_param.is_update_inner_var)
            {
                objID_match_list_new.push_back(_objID_match_list[i]);
                error_list_new.push_back(_error_list[i]);
            }
        }
    }
    if (_param.is_update_inner_var)
    {
        _objID_match_list = objID_match_list_new;
        _error_list = error_list_new;
    }

    // print info
    _n_after_del = tr3d_list.size();
    _n_del = _n_before_del - _n_after_del;
    std::cout << "\tFinish deleting gohst match: "
              << "n_del = " << _n_del << ", "
              << "n_after_del = " << _n_after_del << "."
              << std::endl;
}

// fill tracer info
void StereoMatch::fillTracerInfo (std::vector<Tracer3D>& tr3d_list, std::vector<std::vector<Tracer2D>> const& tr2d_list)
{
    Tracer3D tr3d;
    tr3d._camid_list = _cam_list.useid_list;
    tr3d._n_2d = _n_cam_use;
    tr3d._tr2d_list.resize(_n_cam_use);

    std::vector<Line3D> sight3D_list(_n_cam_use);
    
    int tr2d_id;
    int cam_id;
    for (int i = 0; i < _objID_match_list.size(); i ++)
    {
        for (int id = 0; id < _n_cam_use; id ++)
        {
            tr2d_id = _objID_match_list[i][id];
            cam_id = _cam_list.useid_list[id];

            tr3d._tr2d_list[id]._pt_center = tr2d_list[id][tr2d_id]._pt_center;
            tr3d._tr2d_list[id]._r_px = tr2d_list[id][tr2d_id]._r_px;
            sight3D_list[id] = _cam_list.cam_list[cam_id].lineOfSight(tr3d._tr2d_list[id]._pt_center);
        }

        tr3d._r2d_px = tr3d._tr2d_list[0]._r_px;
        myMATH::triangulation(tr3d._pt_center, tr3d._error, sight3D_list);

        tr3d_list.push_back(tr3d);
    }
}

// save tracer info
void StereoMatch::saveTracerInfo (std::string path, std::vector<Tracer3D> const& tr3d_list)
{
    std::ofstream file;
    file.open(path, std::ios::out);

    file << "WorldX,WorldY,WorldZ,Error,Ncam";

    int n_cam_all = _cam_list.cam_list.size();
    for (int i = 0; i < n_cam_all; i ++)
    {
        file << "," << "cam" << i << "_" << "x(col)" 
             << "," << "cam" << i << "_" << "y(row)";
    }
    file << "\n";

    file.precision(SAVEPRECISION);

    for (int i = 0; i < tr3d_list.size(); i ++)
    {
        tr3d_list[i].saveObject3D(file, n_cam_all);
    }

    file.close();

}

// auxiliary functions for tracerMatch //

// recursively find matches for tracer
void StereoMatch::findTracerMatch (
    int id, // id = 1 => 2nd used cam 
    std::vector<int> const& trID_match,
    std::deque<std::vector<int>>& trID_match_list,
    std::deque<double>& error_list,
    std::vector<std::vector<Tracer2D>> const& tr2d_list
)
{
    if (id < 1)
    {
        std::cerr << "StereoMatch::findTracerMatch error at line " << __LINE__ << ":\n"
                  << "id = " << id << " < 1"
                  << std::endl;
        throw error_size;
    }

    int camID_curr = _cam_list.useid_list[id];
    int n_row = _cam_list.cam_list[camID_curr].getNRow();
    int n_col = _cam_list.cam_list[camID_curr].getNCol();

    // Line of sight
    std::vector<Line2D> sight2D_list;
    std::vector<Line3D> sight3D_list;
    Line2D sight2D;
    Line3D sight3D;
    Pt2D pt2d_1;
    Pt2D pt2d_2;
    Pt2D unit2d;
    for (int i = 0; i < id; i ++)
    {       
        // project from cam_prev onto cam_curr
        int camID_prev = _cam_list.useid_list[i];

        // get 3d light of sight from cam_prev
        sight3D = _cam_list.cam_list[camID_prev].lineOfSight(tr2d_list[i][trID_match[i]]._pt_center);
        sight3D_list.push_back(sight3D);

        // project 3d light of sight onto cam_curr (3d line => 2d line)
        pt2d_1 = _cam_list.cam_list[camID_curr].project(sight3D.pt);
        pt2d_2 = _cam_list.cam_list[camID_curr].project(sight3D.pt + sight3D.unit_vector);
        unit2d = myMATH::createUnitVector(pt2d_1, pt2d_2);
        sight2D.pt = pt2d_1;
        sight2D.unit_vector = unit2d;
        sight2D_list.push_back(sight2D);
    }
    sight3D_list.push_back(sight3D);


    //                       //
    // id = 1 (2nd used cam) //
    //                       //
    Pt3D pt3d;
    if (id == 1)  
    {
        // find tracer candidates around the line of sight         
        // ObjectMarkerAroundLine   

        // calculate the parallel border lines: represented by the ref points
        Line2D sight2D_plus;
        Line2D sight2D_minus;
        // calculate the perpendicular unit vecotr
        unit2d[0] = sight2D.unit_vector[1];
        unit2d[1] = -sight2D.unit_vector[0];
        // calculate the ref points
        pt2d_1 = sight2D.pt + unit2d * _param.tor_2d;
        pt2d_2 = sight2D.pt - unit2d * _param.tor_2d;
        sight2D_plus.pt = pt2d_1;
        sight2D_plus.unit_vector = sight2D.unit_vector;
        sight2D_minus.pt = pt2d_2;
        sight2D_minus.unit_vector = sight2D.unit_vector;

        // Search the candidate within the torlerance near a line
        // Determine the border of the torlerance by iterating each x pixel or y pixel depending on the slope.
        // x_pixel (col id), y_pixel (row id)
        //  if the |slope| > 1, iterate every y_pixel (row)
        //  else, iterate every x_pixel (col) 
        Line2D sight2D_axis;
        if (std::fabs(sight2D.unit_vector[1]) > std::fabs(sight2D.unit_vector[0]))
        {
            sight2D_axis.unit_vector = Pt2D(1,0);
            int x_pixel_1, x_pixel_2, min, max;   

            int y_start = 0;
            int y_end = n_row;
            for (int y_pixel = y_start; y_pixel < y_end; y_pixel ++)
            {
                sight2D_axis.pt[0] = 0;
                sight2D_axis.pt[1] = y_pixel;

                // Get x_pixel (col id) from two lines 
                pt2d_1 = myMATH::crossPoint(sight2D_axis, sight2D_plus);
                pt2d_2 = myMATH::crossPoint(sight2D_axis, sight2D_minus);
                x_pixel_1 = pt2d_1[0];
                x_pixel_2 = pt2d_2[0];
                if (x_pixel_1 > x_pixel_2)
                {
                    max = x_pixel_1 + 1; 
                    min = x_pixel_2;
                }
                else 
                {
                    max = x_pixel_2 + 1; 
                    min = x_pixel_1;
                }
                min = std::max(0, min);
                min = std::min(n_col-1, min);
                max = std::max(1, max);
                max = std::min(n_col, max);

                for (int col = min; col < max; col ++)
                {
                    iterOnObjIDMap (
                        id, 
                        y_pixel, col,
                        sight2D_list,
                        sight3D_list,
                        trID_match, 
                        trID_match_list,
                        error_list,
                        tr2d_list
                    );
                }
            }
        }
        else 
        {
            sight2D_axis.unit_vector = Pt2D(0,1);
            int y_pixel_1, y_pixel_2, min, max;

            int x_start = 0;
            int x_end = n_col;
            for (int x_pixel = x_start; x_pixel < x_end; x_pixel ++)
            {
                sight2D_axis.pt[0] = x_pixel;
                sight2D_axis.pt[1] = 0;

                // Get y_pixel (row id) from two lines 
                pt2d_1 = myMATH::crossPoint(sight2D_axis, sight2D_plus);
                pt2d_2 = myMATH::crossPoint(sight2D_axis, sight2D_minus);
                y_pixel_1 = pt2d_1[1];
                y_pixel_2 = pt2d_2[1];
                if (y_pixel_1 > y_pixel_2)
                {
                    max = y_pixel_1 + 1; 
                    min = y_pixel_2;
                }
                else 
                {
                    max = y_pixel_2 + 1; 
                    min = y_pixel_1;
                }
                min = std::max(0, min);
                min = std::min(n_row-1, min);
                max = std::max(0, max);
                max = std::min(n_row, max);

                for (int row = min; row < max; row ++)
                {
                    iterOnObjIDMap (
                        id, 
                        row, x_pixel,
                        sight2D_list,
                        sight3D_list,
                        trID_match, 
                        trID_match_list,
                        error_list,
                        tr2d_list
                    );
                }
            }
        }

    }
    //                                  //
    // cam_cur_id = 2 ~ _check_cam_id-1 //
    //                                  //
    else
    {
        // find search region
        PixelRange search_region = findSearchRegion(id, sight2D_list);

        // iterate every pixel in the search region
        for (int i = search_region.row_min; i < search_region.row_max; i ++)
        {
            for (int j = search_region.col_min; j < search_region.col_max; j ++)
            {
                // judge whether the distances between the candidate
                // and all the lines are all within the range 
                iterOnObjIDMap (
                    id, 
                    i, j,
                    sight2D_list,
                    sight3D_list,
                    trID_match, 
                    trID_match_list,
                    error_list,
                    tr2d_list
                );
            }
        }
    }
}


void StereoMatch::iterOnObjIDMap (
    int id, 
    int row_id, int col_id,
    std::vector<Line2D> const& sight2D_list,
    std::vector<Line3D>& sight3D_list,
    std::vector<int> const& trID_match, 
    std::deque<std::vector<int>>& trID_match_list,
    std::deque<double>& error_list,
    std::vector<std::vector<Tracer2D>> const& tr2d_list
)
{
    int camID_curr = _cam_list.useid_list[id];
    Pt3D pt3d;
    double tor_2d_sqr = _param.tor_2d * _param.tor_2d;

    for (
        int k = 0; 
        k < _objID_map_list[id](row_id, col_id).size(); 
        k ++
    )
    {
        int tr_id = _objID_map_list[id](row_id, col_id)[k];

        if (tr_id == -1)
        {
            break;
        }

        bool in_range = true;
        for (int m = 0; m < id; m ++)
        {
            double dist2 = myMATH::dist2(tr2d_list[id][tr_id]._pt_center, sight2D_list[m]);
            if (dist2 > tor_2d_sqr)
            {
                in_range = false;
                break;
            }
        } 

        // reproject onto the previous cam, then is within error line
        if (in_range && checkReProject(id, tr_id, trID_match, tr2d_list))
        {
            std::vector<int> trID_match_new = trID_match;
            trID_match_new.push_back(tr_id);

            // if there is still other cameras to check
            if (id < _param.check_id - 1) 
            {
                int next_id = id + 1;

                // move onto the next camera and search candidates
                findTracerMatch (
                    next_id,
                    trID_match_new, 
                    trID_match_list,
                    error_list,
                    tr2d_list                  
                );
            }
            // if the current camera is the last one, then finalize the match.
            else
            {
                // test 3d distance by triangulation  
                sight3D_list[id] = _cam_list.cam_list[camID_curr].lineOfSight(tr2d_list[id][tr_id]._pt_center);
                double error_3d = 0.0;
                myMATH::triangulation(pt3d, error_3d, sight3D_list);

                if (error_3d > _param.tor_3d)
                {
                    continue;
                }
                else 
                {
                    if (id < _n_cam_use - 1)
                    {
                        int next_id = id + 1;

                        checkTracerMatch (
                            next_id,
                            pt3d,
                            trID_match_new,
                            trID_match_list,
                            error_list,
                            tr2d_list
                        );
                    }
                    else 
                    {
                        trID_match_list.push_back(trID_match_new);

                        error_list.push_back(error_3d);
                    }
                }
            }
        }

    }  
}


PixelRange StereoMatch::findSearchRegion (int id, std::vector<Line2D> const& sight2D_list)
{
    PixelRange search_region;

    int n_row = _cam_list.cam_list[_cam_list.useid_list[id]].getNRow();
    int n_col = _cam_list.cam_list[_cam_list.useid_list[id]].getNCol();

    if (sight2D_list.size() == 1)
    {
        search_region.row_max = n_row;
        search_region.col_max = n_col;
        return search_region;
    }

    // Find all crossing points
    Pt2D pt_cross;
    Pt2D perp_unit_curr;
    Pt2D perp_unit_test;
    std::vector<Line2D> line_curr_list(2);
    std::vector<Line2D> line_test_list(2);

    int search_region_state = 0;
    for (int i = 0; i < id-1; i ++)
    {
        perp_unit_curr[0] = sight2D_list[i].unit_vector[1];
        perp_unit_curr[1] = - sight2D_list[i].unit_vector[0];

        for (int j = i+1; j < id; j ++)
        {
            perp_unit_test[0] = sight2D_list[j].unit_vector[1];
            perp_unit_test[1] = - sight2D_list[j].unit_vector[0];

            // find the bounded lines for each crossing line
            line_curr_list[0].pt = sight2D_list[i].pt + perp_unit_curr * _param.tor_2d;
            line_curr_list[0].unit_vector = sight2D_list[i].unit_vector;
            line_curr_list[1].pt = sight2D_list[i].pt - perp_unit_curr * _param.tor_2d;
            line_curr_list[1].unit_vector = sight2D_list[i].unit_vector;

            line_test_list[0].pt = sight2D_list[j].pt + perp_unit_test * _param.tor_2d;
            line_test_list[0].unit_vector = sight2D_list[j].unit_vector;
            line_test_list[1].pt = sight2D_list[j].pt - perp_unit_test * _param.tor_2d;
            line_test_list[1].unit_vector = sight2D_list[j].unit_vector;

            // find the crossing points
            for (int k = 0; k < 2; k ++)
            {
                for (int l = 0; l < 2; l ++)
                {
                    pt_cross = myMATH::crossPoint(line_curr_list[k], line_test_list[l]);

                    if (search_region_state == 0)
                    {
                        search_region.row_min = int(pt_cross[1]);
                        search_region.row_max = search_region.row_min + 1;
                        search_region.col_min = int(pt_cross[0]);
                        search_region.col_max = search_region.col_min + 1;

                        search_region_state = 1;
                    }
                    else 
                    {
                        search_region.setRange(pt_cross[1], pt_cross[0]);
                    }
                }
            }
        }
    }

    search_region.row_min = std::min(std::max(search_region.row_min, 0), n_row-1);
    search_region.row_max = std::min(std::max(search_region.row_max, 1), n_row);
    search_region.col_min = std::min(std::max(search_region.col_min, 0), n_col-1);
    search_region.col_max = std::min(std::max(search_region.col_max, 1), n_col);

    return search_region;
}


bool StereoMatch::checkReProject (
    int id, 
    int tr_id,
    std::vector<int> const& tracer_id_match, 
    std::vector<std::vector<Tracer2D>> const& tr2d_list
)
{
    int camID_curr = _cam_list.useid_list[id];
    Line3D sight3D = _cam_list.cam_list[camID_curr].lineOfSight(tr2d_list[id][tr_id]._pt_center);

    Line2D sight2D;
    Pt2D pt2d_1;
    Pt2D pt2d_2;
    double dist2 = 0;
    double tor_2d_sqr = _param.tor_2d * _param.tor_2d;

    for (int i = 0; i < id; i ++)
    {
        int cam_id = _cam_list.useid_list[i];

        pt2d_1 = _cam_list.cam_list[cam_id].project(sight3D.pt);
        pt2d_2 = _cam_list.cam_list[cam_id].project(sight3D.pt + sight3D.unit_vector);
        sight2D.pt = pt2d_1;
        sight2D.unit_vector = myMATH::createUnitVector(pt2d_1, pt2d_2);

        dist2 = myMATH::dist2(tr2d_list[i][tracer_id_match[i]]._pt_center, sight2D);

        if (dist2 > tor_2d_sqr)
        {
            return false;
        }
    }

    return true;
}


void StereoMatch::checkTracerMatch(
    int id, 
    Pt3D const& pt3d,
    std::vector<int> const& trID_match, 
    std::deque<std::vector<int>>& trID_match_list,
    std::deque<double>& error_list,
    std::vector<std::vector<Tracer2D>> const& tr2d_list
)
{
    if (id < 1)
    {
        std::cerr << "StereoMatch::checkTracerMatch error at line " << __LINE__ << ":\n"
                  << "id = " << id << " < 1"
                  << std::endl;
        throw error_size;
    }

    int camID_curr = _cam_list.useid_list[id];
    int n_row = _cam_list.cam_list[camID_curr].getNRow();
    int n_col = _cam_list.cam_list[camID_curr].getNCol();

    // Line of sight
    std::vector<Line2D> sight2D_list;
    std::vector<Line3D> sight3D_list;
    Line2D sight2D;
    Line3D sight3D;
    Pt2D pt2d_1;
    Pt2D pt2d_2;
    Pt2D unit2d;
    for (int i = 0; i < id; i ++)
    {       
        // project from cam_prev onto cam_curr
        int camID_prev = _cam_list.useid_list[i];

        // get 3d light of sight from cam_prev
        sight3D = _cam_list.cam_list[camID_prev].lineOfSight(tr2d_list[i][trID_match[i]]._pt_center);
        sight3D_list.push_back(sight3D);

        // project 3d light of sight onto cam_curr (3d line => 2d line)
        pt2d_1 = _cam_list.cam_list[camID_curr].project(sight3D.pt);
        pt2d_2 = _cam_list.cam_list[camID_curr].project(sight3D.pt + sight3D.unit_vector);
        unit2d = myMATH::createUnitVector(pt2d_1, pt2d_2);
        sight2D.pt = pt2d_1;
        sight2D.unit_vector = unit2d;
        sight2D_list.push_back(sight2D);
    }
    sight3D_list.push_back(sight3D);

    // find search region
    //  directly project pt_world onto the image plane of the current camera
    //  then find the search region
    Pt2D pt2d = _cam_list.cam_list[camID_curr].project(pt3d);

    int row_min = pt2d[1] - _param.check_radius;
    int row_max = pt2d[1] + _param.check_radius + 1;
    int col_min = pt2d[0] - _param.check_radius;
    int col_max = pt2d[0] + _param.check_radius + 1;

    row_min = std::max(0, row_min);
    row_max = std::min(n_row, row_max);
    col_min = std::max(0, col_min);
    col_max = std::min(n_col, col_max);

    // if the search region is out of the image plane, then return
    if (row_min >= row_max || col_min >= col_max)
    {
        return;
    }

    for (int i = row_min; i < row_max; i ++)
    {
        for (int j = col_min; j < col_max; j ++)
        {
            iterOnObjIDMap (
                id, 
                i, j,
                sight2D_list,
                sight3D_list,
                trID_match, 
                trID_match_list,
                error_list,
                tr2d_list
            );
        }
    }
}

#endif