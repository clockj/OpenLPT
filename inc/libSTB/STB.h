#ifndef STB_H
#define STB_H

#include "myMATH.h"
#include "Matrix_bool.h"
#include "Matrix.h"
#include "myIO.h"
#include "Camera.h"
#include "ObjectInfo.h"
#include "ImageIO.h"
#include "ObjectFinder.h"
#include "StereoMatch.h"
#include "OTF.h"
#include "Shake.h"
#include "IPR.h"
#include "PredField.h"
#include "Track.h"

#include <vector>
#include <deque>
#include <typeinfo>
#include <omp.h>
#include <ctime>
#include <iostream>
#include <fstream>
#include <sstream>

template<class T>
class STB
{
private:
    static const int _UNLINKED = -1;
    static const int _MINTRACK = 10;
    static const int _EMPTY = INT_MAX;


    // ####### IO data ####### //
    // first and last frames
    int _first, _last;

    // view area
    AxisLimit _xyz_limit;
    double _vox_to_mm;

    // camera param
    int _n_cam, _n_pix_w, _n_pix_h;
    std::string _img_folder;
    std::vector<Camera> _cam_list;
    std::vector<ImageIO> _imgio_list;
    
    // initialization phase
    std::vector<std::vector<T>> _ipr_matched;
    bool _ipr_flag;	        // use ipr (1) or .mat files (0) for 3D positions?
    double _r_search_init;  // mm 

    // convergence phase
    std::vector<Matrix<double>> _img_org_list;			// original images
    // std::deque<Matrix<double>> _img_res_list;			// residue images
    // std::deque<Matrix<double>> _img_prj_list;           // reproject images 
    double _shake_shift;                                // Shaking range, mm
    double _space_avg;                                  // avg interparticle dist
    double _shift_betframe_max;                         // Largest expected particle shift
    double _shift_abs_max;                              // maximum allowed change in particle shift accross 2 frames
    double _shift_rel_max;                              // maximum allowed relative achange in particle shift accross 2 frames
    double _fpt;										// a multiplying factor on particle intensity in order to ensure the residual has no traces of tracked particles
    double _int_lower;                                  // lower intensity threshold (xx*avg. Intensity) to eliminate ghost particles while tracking

    // Predictive Field 
    std::vector<int> _ngrid_xyz;
    double _r_search_pred; // mm

    // IPR
    bool _TRIONLY = false;                     
    bool _IPRONLY = false;
    int _ipr_loop_outer;                       // no. of outerloop iterations
    int _ipr_loop_inner;                       // no. of innerloop iterations
    double _ipr_int_lower;                     // lower intensity threshold (xx*avg. Intensity) to eliminate ghost particles
    double _ipr_tol_2d;                        // mm  
    double _ipr_tol_3d;                        // mm
    bool _ipr_is_reduced = false;
    int _ipr_loop_outer_reduced;
    double _ipr_tol_2d_reduced;
    double _ipr_tol_3d_reduced;
    std::vector<int> _cam_id_list;
    std::deque<Camera> _cam_reduced;
    OTF _otf;

    // Object info
    T _obj_info_sample;                              // object searching particle radius, pixel
    std::vector<int> _imgint_max_list;
    std::vector<int> _imgint_min_list;


    // ####### STB intermediate variable ####### //
    double _acc_diff_thred = 0;
    bool _is_acc_check = false;

    // storing the tracks
    std::deque<Track<T>> _active_long_track;		    // tracks with more than 3 particles
    std::deque<Track<T>> _active_short_track;		    // tracks with 3 or less particles
    std::deque<Track<T>> _inactive_track;				// tracks that are inactive as they could not find the correct link
    std::deque<Track<T>> _exit_track;					// tracks that left the measurement domain (out of at least 2 cameras)
    std::deque<Track<T>> _inactive_long_tracks;
    std::deque<Track<T>> _buffer_tracks;			    //  new long tracks added from Back or Forward STB (multi-pass)

    // dummy variables to identify the no. of tracks added and subtracted 
    int _a_as = 0, _a_al = 0, _a_is = 0, _s_as1 = 0, _s_as2 = 0, _s_as3 = 0, _s_as4 = 0, _s_al = 0, _a_il = 0;
    //TESTING
    std::vector<Matrix<double>> _temp_pred;
    double _r_design = 5;
    double _lfit_error[100];

    void SetTracerParam (std::stringstream& parsed);

public:
    // constructor
    STB(std::string file);
    // // copy constructor
    // STB(STB& s);
    // destructor
    ~STB() {};

    enum TrackType{ Inactive = 0, ActiveShort = 1, ActiveLong = 2, Exit = 3, InactiveLong = 4, Buffer = 5};

    //############################### FUNCTIONS ##############################
    void Run();
    
    // to make tracks for the first four frames
    void InitialPhase();
    void LoadTracks(std::string path, TrackType trackType);

//     // a function to start a track for particles that were left untracked in current frame
    void StartTrack(int frame, PredField<T>& pf);
    // extends the particle track to nextFrame using search radius method (can be used for first 4 links in both, intialization and convergence phase)
    void MakeLink(int nextframe, const Matrix<double>& vel_curr, double radius, Track<T>& track, bool& active);
    // obj_id = -1 => UNLINKED
    void NearestNeighbor(std::vector<T>& obj_list, double radius, const Matrix<double>& pt_estimate, int& obj_id);
    // obj_id = -1 => UNLINKED
    void NearestNeighbor(std::vector<T>& obj_list, double radius, const Matrix<double>& pt_estimate, int& obj_id, std::vector<int>& candidate_used);

    // convergence phase
    void ConvergencePhase();
    void Prediction(int frame, std::vector<Matrix<double>>& est_pos);
//     vector<double> Polyfit(Track tracks, string direction, int datapoints, int polydegree);		// predictor for convergence phase
//     /*
//      * Function: predict the next point with Wiener Predictor using LMS algorithm
//      * Input: tracks: the track to be predicted
//      * 		  direction: to indicate the axis to be predicted
//      * 		  order: the order of the Wiener Predictor
//      * Notice: the order should be less than the length of track.
//      * Output: the next point
//      */
    double LMSWienerPred(Track<T>& track, std::string direction, int order);

    // Link short tracks with obj candidates in residual images
    void MakeShortLinkResidual(int nextframe, std::vector<T>& obj_list, Track<T>& tr, int n_iter, int& is_erase, std::vector<int>& candidate_used);

// //	bool CheckVelocity(deque<Track>::iterator& tr);
//     bool CheckVelocity(Track& tr);
//     bool CheckAcceleration(deque<Track>::iterator& tr);
//     void GetAccThred();
    bool CheckLinearFit(Track<T>& tr);

//     void MatTracksSave(string addres, string s, bool is_back_STB);
//     void MatfileSave(deque<Track> tracks, string address, string name, int size);
//     void SaveTrackToTXT(deque<Track> tracks, string address);
//     void LoadAllTracks(string address, string frame_number, bool is_back_STB);
//     void LoadTrackFromTXT(string path, TrackType trackType);

//     //###################### TEMPORARY FUNTIONS FOR TESTING ###############################

//     // to load 3D positions without IPR
//     std::vector<Matrix<double>> Load_3Dpoints(string path);

};

#include "STB.hpp"

#endif