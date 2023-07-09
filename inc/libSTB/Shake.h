#ifndef SHAKE_H
#define SHAKE_H

#include <vector>
#include <string>
#include <typeinfo>
#include <omp.h>

#include "Matrix.h"
#include "ObjectInfo.h"
#include "Camera.h"
#include "OTF.h"


template<class T>
class Shake
{
protected:
    // INPUTS //
    std::vector<T>& _object_info;               // 3D
    std::vector<Matrix<double>>& _orig_img_list;   // Original images
    double _shake_width;                        // unit: mm
    std::vector<Camera>& _cam;                  // Camera parameter
    int _n_cam;                                 // Number of camera
    OTF& _otf;                                  // Camera OTF
    std::vector<int> _max_intensity;            // default: 8 digit (255)
    int _n_loop = 4;                            // Number of shake times
    double _int_low_ratio = 0.1;                // Ghost threshold

    std::vector<int> _cam_id;                   // Cam id that is used for shaking (so that no need to rearrange _orig_img_list)
    std::vector<Camera> _cam_use;               // Cam that is in use
    int _n_cam_use;                             // Number of cam in use


    // INTERMEDIATES //
    std::vector<Matrix<double>> _res_img_list; // Original image - Particle projection image
    std::vector<double> _intensity_list;


    // OUTPUTS //
    std::vector<int> _remove_object_id;
    int _n_ghost = 0;


    // FUNCTIONS //
    // Calculate the windows size for shaking
    PixelRange FindSearchRegion (int cam_use_id, int center_row, int center_col, int half_width_pixel); // a square region
    
    // Projection for particles
    double GaussianProjection (Matrix<double>& tracer_center_pixel, std::vector<double>& otf_param, int x, int y); // (x,y)=(col,row)

    // Calculate residue for shaking
    double CalPointResidue (Matrix<double>& pos_new, std::vector<PixelRange>& search_range_list, std::vector<Matrix<double>>& aug_img_list);

    // Search for the minimum residue
    void FitQuadratic (std::vector<double>& coeff, double* array, std::vector<double>& residue);

    // Shaking and refine 3D position
    Matrix<double> UpdatePos3D (Matrix<double> const& pos_old, std::vector<PixelRange>& search_range_list, std::vector<Matrix<double>>& aug_img_list, double delta);

    // Calculate intensity for shaken particles
    double CalObjectIntensity (TracerInfo& tracer, std::vector<PixelRange>& search_range_list, std::vector<Matrix<double>>& aug_img_list, double intensity, double min_ratio_of_median);

    // Procedure for each shake
    double ShakingTracer(TracerInfo& tracer, double delta, double intensity);

    // Remove all tracked particles from image to get residual image.
    void TracerResImg();

    // Remove negative pxiel and set them as zeros, this function is used to prepare residual image for the next run of IPR.
    void AbsResImg ();

    // Remove ghost particles.
    void RemoveGhost();

public:
// pending: 07/05 use cam_id
    Shake(
        std::vector<Matrix<double>>& orig_img_list, // all img list
        std::vector<T>& object_info, // also used as output
        double shake_width, // unit: mm
        std::vector<Camera>& cam, // all cam parameter
        OTF& otf // all cam otf parameter
    ) : _orig_img_list(orig_img_list), _object_info(object_info), _shake_width(shake_width), _cam(cam), _n_cam(cam.size()), _otf(otf), _res_img_list(orig_img_list), _intensity_list(object_info.size(), 1)
    {
        // default for max_intensity
        SetMaxIntensity(255);

        // default for cam_id
        _n_cam_use = _n_cam;
        for (int i = 0; i < _n_cam_use; i ++)
        {
            _cam_id.push_back(i);
            _cam_use.push_back(_cam[i]);
        }
    };

    Shake(
        std::vector<Matrix<double>>& orig_img_list, // all img list
        std::vector<T>& object_info, // also used as output
        double shake_width, // unit: mm
        std::vector<Camera>& cam, // all cam parameter
        OTF& otf, // all cam otf parameter
        std::vector<int>& cam_id
    ) : _orig_img_list(orig_img_list), _object_info(object_info), _shake_width(shake_width), _cam(cam), _n_cam(cam.size()), _otf(otf), _cam_id(cam_id), _n_cam_use(cam_id.size()), _res_img_list(orig_img_list), _intensity_list(object_info.size(), 1)
    {
        // default for max_intensity
        SetMaxIntensity(255);

        int id;
        for (int i = 0; i < _n_cam_use; i ++)
        {
            id = _cam_id[i];
            _cam_use.push_back(_cam[id]);
        }
    };

    ~Shake() {};

    void RunShake();

    // For output: Matrix<int> residue img 
    void GetResImg(std::vector<Matrix<double>>& res_img_list_out) 
    {
        int n_cam = _cam.size();
        for (int i = 0; i < n_cam; i ++) 
        {
            res_img_list_out[i] = _res_img_list[i];
        }
    };

    void SetMaxIntensity(int max_intensity)
    {
        _max_intensity = std::vector<int>(_n_cam, max_intensity);
    };

    void SetMaxIntensity(std::vector<int> max_intensity)
    {
        if (max_intensity.size() != _n_cam)
        {
            std::cerr << "The size of input maximum intensity vector is "
                << max_intensity.size() << "; "
                << "it does not match the size of camera" << _cam.size()
                << std::endl;
        }
        _max_intensity = max_intensity;
    };

    void SetShakeTimes (int n_loop) 
    {
        _n_loop = n_loop;
    };

    int GetNGhost () 
    {
        return _n_ghost;
    };

    double GetIntLowRatio ()
    {
        return _int_low_ratio;
    };
    void SetIntLowRatio (double ratio)
    {
        _int_low_ratio = ratio;
    };

    void SetCamID (std::vector<int> cam_id)
    {
        _cam_id = cam_id;
        _n_cam_use = cam_id.size();

        _cam_use = std::vector<Camera> (_n_cam_use);
        for (int m = 0; m < _n_cam_use; m ++)
        {
            _cam_use[m] = _cam[_cam_id[m]];
        }
    };

};

#include "Shake.hpp"

#endif // !SHAKE_H

