//
//	StereoMatch.h 
//
//	Base class for stereomatch
//
//	Created by Shijie Zhong 07/06/2022 
//

#ifndef STEREOMATCH_H
#define STEREOMATCH_H

#include <vector>
#include <omp.h>
#include <ctime>
#include <string>
#include <deque>
#include <algorithm>

#include <iostream>
#include <string>
#include <sstream>
#include <fstream>

#include "ObjectInfo.h"
#include "Camera.h"
#include "Matrix.h"
#include "STBCommons.h"
#include "myMATH.h"

// Map from (row_id,col_id) to object ID for each camera
class ObjIDMap
{
private:
    int _n_row = 0;
    int _n_col = 0;

    // img_id = row_id * _n_col + col_id
    // map from img_id to object ID
    // _map[img_id] = object ID list on (row_id, col_id)
    std::vector<std::vector<int>> _map; 

public:
    ObjIDMap() {};
    ~ObjIDMap() {};

    void config(int n_row, int n_col)
    {
        _n_row = n_row; 
        _n_col = n_col;

        _map.reserve(_n_row * _n_col);
        _map.resize(_n_row * _n_col);
        std::fill(_map.begin(), _map.end(), std::vector<int>(1,-1));
    };

    // map from (row_id, col_id) to img_id
    int mapImgID (int row_id, int col_id) 
    {
        return row_id * _n_col + col_id;
    };

    std::vector<int>& operator() (int row_id, int col_id) 
    {
        return _map[mapImgID(row_id, col_id)];
    };
};


// Stereo match parameters
struct SMParam
{
    double tor_2d = -1;   // [px] 2D tolerance
    double tor_3d = -1;   // [mm] 3D tolerance
    int n_thread = 0;     // num of parallel threads
    int check_id = 2;     // check_id <= n_use !
    double check_radius = 1; // [px]
    bool is_delete_ghost = false; // delete ghost tracer
    bool is_update_inner_var = false; // update inner variables
};


class StereoMatch
{
public:
    /**************************INPUT VARIABLES******************************/
    SMParam _param;
    CamList const& _cam_list; // camera list
    int _n_cam_use = 0; // number of cameras used

    /*************************OUTPUT VARIABLES******************************/
    std::vector<std::vector<int>> _objID_match_list; // matches of object ID
    std::vector<double> _error_list; // error list

    StereoMatch(SMParam const& param, CamList const& cam_list);
    ~StereoMatch() {};

    void clearAll ();

    /*******************************MAIN**********************************/
    // Match 3D object to 2D object
    // obj3d_list: 3D object list (output)
    // obj2d_list: 2D object list (input)
    template<class T3D, class T2D>
    void match(std::vector<T3D>& obj3d_list, std::vector<std::vector<T2D>> const& obj2d_list);

    template<class T3D>
    void saveObjInfo (std::string path, std::vector<T3D> const& obj3d_list);
    void saveObjIDMatchList (std::string path);
    
private:
    /*************************PROCESS VARIABLES*****************************/
    std::vector<ObjIDMap> _objID_map_list;  // image map for objects

    int _n_del = 0;
    int _n_before_del = 0;
    int _n_after_del = 0;


    /*****************************Functions*********************************/
    // Create object ID map
    // This map facilitates the search of object on a specific image position
    template<class T2D>
    void createObjIDMap (std::vector<std::vector<T2D>> const& obj2d_list);

    //         //
    // Tracers //
    //         //
    void tracerMatch(std::vector<std::vector<Tracer2D>> const& tr2d_list);
    
    void removeGhostTracer (std::vector<Tracer3D>& obj3d_list, std::vector<std::vector<Tracer2D>> const& tr2d_list);
    void removeGhostTracerTest (std::vector<Tracer3D>& obj3d_list, std::vector<std::vector<Tracer2D>> const& tr2d_list);

    void fillTracerInfo (std::vector<Tracer3D>& tr3d_list, std::vector<std::vector<Tracer2D>> const& tr2d_list);
    
    void saveTracerInfo (std::string path, std::vector<Tracer3D> const& tr3d_list);

    // recursively find matches for tracer
    // trID_match_list: all matches for tracer A in camera 1 after going through all cameras
    // trID_match: current match for tracer A till the current camera
    void findTracerMatch (
        int id,
        std::vector<int> const& trID_match,
        std::deque<std::vector<int>>& trID_match_list, 
        std::deque<double>& error_list,
        std::vector<std::vector<Tracer2D>> const& tr2d_list
    );

    PixelRange findSearchRegion (int id, std::vector<Line2D> const& sight2D_list);

    void iterOnObjIDMap (
        int id, 
        int row_id, int col_id,
        std::vector<Line2D> const& sight2D_list,
        std::vector<Line3D>& sight3D_list,
        std::vector<int> const& trID_match, 
        std::deque<std::vector<int>>& trID_match_list,
        std::deque<double>& error_list,
        std::vector<std::vector<Tracer2D>> const& tr2d_list
    );

    bool checkReProject(
        int id, 
        int tr_id,
        std::vector<int> const& trID_match, 
        std::vector<std::vector<Tracer2D>> const& tr2d_list
    );

    void checkTracerMatch(
        int id, 
        Pt3D const& pt3d,
        std::vector<int> const& trID_match,
        std::deque<std::vector<int>>& trID_match_list,
        std::deque<double>& error_list,
        std::vector<std::vector<Tracer2D>> const& tr2d_list
    );

};

#include "StereoMatch.hpp"

#endif

