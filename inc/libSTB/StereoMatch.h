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
struct StereoMatchParam
{
    double search_radius = 0; // [px] search radius
    double tor_2d = -1; // [px] 2D tolerance
    double tor_3d = -1; // [mm] 3D tolerance
    int n_thread = 0;   // num of parallel threads
    int check_id; // _check_cam_id <= n_use
    double fov_size = 40; // [mm] 
    bool is_delete_ghost = false; // delete ghost tracer
};


class StereoMatch
{
public:
    /**************************INPUT VARIABLES******************************/
    StereoMatchParam _param;
    CamList& _cam_list; // camera list
    int _n_cam_use = 0; // number of cameras used

    /*************************OUTPUT VARIABLES******************************/
    std::vector<double> _error_list; // error list

    StereoMatch(StereoMatchParam param, CamList& cam_list);
    ~StereoMatch() {};

    void writeTracer2DList (std::string file, std::vector<std::vector<Tracer2D>> const& tr2d_list);
    void writeTriErrorList (std::string file);
    void writeObjMatchList (std::string file);

    void clearAll ();

    /*******************************MAIN**********************************/
    // Create object ID map
    // This map facilitates the search of object on a specific image position
    template<class T>
    void createObjIDMap (std::vector<std::vector<T>> const& obj2d_list);

    // Match 3D object to 2D object
    // obj3d_list: 3D object list (output)
    // obj2d_list: 2D object list (input)
    template<class T>
    void match(std::vector<T>& obj3d_list, std::vector<std::vector<T>> const& obj2d_list);
    
private:
    /*************************PROCESS VARIABLES*****************************/
    std::vector<ObjIDMap> _objID_map_list;  // image map for objects
    std::vector<std::vector<int>> _objID_match_list; // matches of object ID

    int _n_del = 0;
    int _n_before_del = 0;
    int _n_after_del = 0;

    /*****************************Functions*********************************/
    template<class T>
    void fillObjectInfo (std::vector<std::vector<T>>& object_list_pixel);


    //         //
    // Tracers //
    //         //
    void tracerMatch(std::vector<Tracer3D>& tr3d_list, std::vector<std::vector<Tracer2D>> const& tr2d_list);

    // Find match for tracer
    // tracer_id_match_list: all matches for tracer A in camera 1 after going through all cameras
    // tracer_id_match: current match for tracer A till the current camera
    void findTracerMatch (
        int id,
        std::vector<int> trID_match,
        std::deque<std::vector<int>>& trID_match_list, 
        std::deque<double>& error_list,
        std::vector<std::vector<Tracer2D>> const& tr2d_list
    );

    PixelRange findSearchRegionOrig (int id, std::vector<Line2D> const& sight2D_list);
    PixelRange findSearchRegion (int id, std::vector<Line2D> const& sight2D_list);
    PixelRange findSearchRegion (int id, int line_ref_id, std::vector<Matrix<double>>& pt_list, std::vector<Matrix<double>>& unit_list);

    void iterOnObjIDMap (
        int id, 
        int row_id, int col_id,
        std::vector<Line2D> const& sight2D_list,
        std::vector<Line3D>& sight3D_list,
        std::vector<int>& trID_match, 
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
        std::deque<std::vector<int>>& trID_match_list,
        std::vector<int>& trID_match, 
        Pt3D& pt3d,
        std::deque<double>& error_list
    );

    int isWithinRange (int cam_cur_id, Matrix<double>& object_pos, std::vector<Matrix<double>>& pt_list, std::vector<Matrix<double>>& line_of_sight);

    void deleteGohstMatch (std::vector<std::vector<Tracer2D>>& tr2d_list);
    void deleteGohstMatchNew (std::vector<std::vector<Tracer2D>>& tr2d_list);
    
    void merge (std::vector<int>& A, int id_begin, int id_middle, int id_end, std::vector<int>& B, std::vector<int>& C);
    void splitMerge (std::vector<int>& B, int id_begin, int i_end, std::vector<int>& A, std::vector<int>& C);
    void sortPreference (std::vector<int>& sorted_id_list, std::vector<int>& object_id_match_rank);
    void deleteGohstMatchOrig (std::vector<std::vector<Tracer2D>>& object_list_pixel);

    
 
};

#include "StereoMatch.hpp"

#endif

