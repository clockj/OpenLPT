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

#include <iostream>
#include <string>
#include <sstream>
#include <fstream>

#include "ObjectInfo.h"
#include "Camera.h"
#include "Matrix.h"
#include "STBCommons.h"
#include "myMATH.h"


// template<ObjectType T>
template<class T>
class StereoMatch
{
protected:
    /**************************INPUT VARIABLES******************************/
    std::vector<Camera> _cam_list;   // list of camera parameters
    int _n_cam = 0;
    double _tor_2d = -1;
    double _tor_3d = -1;
    int _check_cam_id; // _check_cam_id <= n_cam
    int _tor_1d = 3;
    
    int _n_thread = 0; // num of parallel threads
    double _fov_length_mm = 40; 

    /*************************PROCESS VARIABLES*****************************/
    std::vector<std::vector<T>> _object_list_mm;  // list of 2D image position in mm
    std::vector<
        std::vector<
            std::vector<
                std::vector<int>
            >
        >
    > _object_id_map;  // image map for objects
    std::vector<std::vector<int>> _object_id_match_list; // matches of object ID
    std::vector<double> _error_list;

    /*************************OUTPUT VARIABLES******************************/
    std::vector<T> _object_info_match_list;
    int _n_del = 0;
    int _n_before_del = 0;
    int _n_after_del = 0;
// public:
    /*****************************Functions*********************************/
    Matrix<double> FindCrossPoint (Matrix<double> const& pt_1, Matrix<double> const& line_of_sight_1, Matrix<double> const& pt_2, Matrix<double> const& line_of_sight_2);
    int ComputeXpixel (Camera& cam, Matrix<double>& pt, Matrix<double>& line_of_sight, double y_pixel);
    int ComputeYpixel (Camera& cam, Matrix<double>& pt, Matrix<double>& line_of_sight, double x_pixel);

    PixelRange FindSearchRegionOrig (int cam_id, std::vector<Matrix<double>>& pt_list, std::vector<Matrix<double>>& unit_list);
    PixelRange FindSearchRegion (int cam_id, std::vector<Matrix<double>>& pt_list, std::vector<Matrix<double>>& unit_list);
    PixelRange FindSearchRegion (int cam_id, int line_ref_id, std::vector<Matrix<double>>& pt_list, std::vector<Matrix<double>>& unit_list);

    void Triangulation(int n_cam_tri, std::vector<Matrix<double>>& pt_2d_list, Matrix<double>& pt_world, double& error);

    // Mark tracers that are within the projection of line of sight.
    // (pt + line of sight) describe the line of sight
    void ObjectMarkerAroundLine (int cam_id, Matrix<double>& pt, Matrix<double>& line_of_sight, Matrix<int>& object_marker, PixelRange& search_region);

    int IsWithinRange (int cam_cur_id, Matrix<double>& object_pos, std::vector<Matrix<double>>& pt_list, std::vector<Matrix<double>>& line_of_sight);

    void FillObjectInfo (std::vector<std::vector<T>>& object_list_pixel);

public:
    StereoMatch(
        std::vector<Camera> const& cam_list,
        double tor_2d,
        double tor_3d
    );

    ~StereoMatch() {};

    /*****************************SET & GET********************************/
    void SetNumThread (int n_thread) {_n_thread = n_thread;};
    void SetTolerance1D (int tor_1d) {_tor_1d = tor_1d;};
    void SetTolerance2D (double tor_2d) {_tor_2d = tor_2d;};
    void SetTolerance3D (double tor_3d) {_tor_3d = tor_3d;};
    void SetCamList (std::vector<Camera> cam_list) {_cam_list = cam_list; _n_cam = cam_list.size();};
    void SetObjectListMMFromMM(std::vector<std::vector<T>>& object_list_mm);

    int GetNumThread () {return _n_thread;};
    int GetTolerance1D () {return _tor_1d;};
    double GetTolerance2D () {return _tor_2d;};
    double GetTolerance3D () {return _tor_3d;};
    std::vector<std::vector<T>> GetObjectListMM () {return _object_list_mm;};
    
    void GetObjectIDMatchList (std::vector<std::vector<int>>& match_list)
    {
        match_list = std::move(_object_id_match_list);
    };
    void GetObjectInfoMatchList (std::vector<T>& object_info_match_list)
    {
        object_info_match_list = std::move(_object_info_match_list);
    };
    void GetTriErrorList (std::vector<double>& error_list)
    {
        error_list = std::move(_error_list);
    };

    void WriteMatchList (std::string file, std::vector<std::vector<T>> const& object_list_pixel);
    void WriteTriErrorList (std::string file);
    void WriteObjectInfoMatchList (std::string file);

    int GetNumOfDelMatch  () {return _n_del;};
    int GetNumOfBeforeDel () {return _n_before_del;};
    int GetNumOfAfterDel  () {return _n_after_del;};

    void ClearCamList      ()      {_cam_list.clear(); _n_cam = 0;};
    void ClearObjectListMM ()      {_object_list_mm.clear();};
    void ClearObjectIDMap  ()      {_object_id_map.clear();};
    void ClearObjectIDMatchList () {_object_id_match_list.clear();};
    void ClearAll ();

    /*******************************MAIN**********************************/
    // Convert the unit of Pixel to MM for the 2D image positions
    void SetObjectListMMFromPixel(std::vector<std::vector<T>>& object_list_pixel);

    // Map object onto an image vector
    // This map facilitates the search of object on a specific image position
    void MakeObjectIDMap (std::vector<std::vector<T>>& object_list_pixel);

    //  input: max_dist_2D (threshold distance between a line of sight and real particles on each image plane [mm])
    void Match(std::vector<std::vector<T>>& object_list_pixel, int is_delete_ghost);
    
private:
    // Find match for tracer
    // tracer_id_match_list: all matches for tracer A in camera 1 after going through all cameras
    // tracer_id_match: current match for tracer A till the current camera
    void FindTracerMatch (
        int cam_cur_id, 
        std::deque<std::vector<int>>& tracer_id_match_list, 
        std::vector<int> tracer_id_match, 
        std::deque<double>& error_list
    );
    void FindTracerMatchOrig (
        int cam_cur_id,
        std::deque<std::vector<int>>& tracer_id_match_list, 
        std::vector<int> tracer_id_match,
        std::deque<double>& error_list
    );
    void CheckTracerMatch(
        int cam_cur_id, 
        std::deque<std::vector<int>>& tracer_id_match_list,
        std::vector<int>& tracer_id_match, 
        Matrix<double>& pt_world,
        std::deque<double>& error_list
    );
    bool CheckReProject(
        int cam_cur_id, 
        std::vector<int> const& tracer_id_match, 
        Matrix<double> const& pt_mm
    );

    void DeleteGohstTracerMatch (std::vector<std::vector<TracerInfo>>& object_list_pixel);
    void DeleteGohstTracerMatchNew (std::vector<std::vector<TracerInfo>>& object_list_pixel);
    
    void Merge (std::vector<int>& A, int id_begin, int id_middle, int id_end, std::vector<int>& B, std::vector<int>& C);
    void SplitMerge (std::vector<int>& B, int id_begin, int i_end, std::vector<int>& A, std::vector<int>& C);
    void SortPreference (std::vector<int>& sorted_id_list, std::vector<int>& object_id_match_rank);
    void DeleteGohstTracerMatchOrig (std::vector<std::vector<TracerInfo>>& object_list_pixel);

    void TracerMatch(std::vector<std::vector<TracerInfo>>& object_list_pixel);
 
};

#include "StereoMatch.hpp"

#endif

