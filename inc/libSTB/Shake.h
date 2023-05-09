#ifndef SHAKE_H
#define SHAKE_H

#include <vector>
#include <string>
#include <typeinfo>

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
    std::vector<Matrix<int>>& _orig_img_list;   // Original images
    double _shake_width;                        // unit: mm
    std::vector<Camera>& _cam;                  // Camera parameter
    OTF& _otf;                                  // Camera OTF
    std::vector<int> _max_intensity;            // default: 8 digit (255)
    int _n_loop = 6;                            // Number of shake times
    double _int_low_ratio = 0.1;                // Ghost threshold


    // INTERMEDIATES //
    std::vector<Matrix<double>> _res_img_list; // Original image - Particle projection image
    std::vector<double> _intensity_list;


    // OUTPUTS //
    std::vector<int> _remove_object_id;
    int _n_ghost = 0;


    // FUNCTIONS //
    // Calculate the windows size for shaking
    PixelRange FindSearchRegion (int cam_id, int center_row, int center_col, int half_width_pixel); // a square region
    
    // Projection for particles
    double GaussianProjection (Matrix<double>& tracer_center_pixel, std::vector<double>& otf_param, int x, int y); // (x,y)=(col,row)

    // Calculate residue for shaking
    double CalPointResidue (Matrix<double>& pos_new, std::vector<PixelRange>& search_range_list, std::vector<Matrix<double>>& aug_img_list);

    // Search for the minimum residue
    void FitQuadratic (std::vector<double>& coeff, double* array, std::vector<double>& residue);

    // Shaking and refine 3D position
    Matrix<double> UpdatePos3D (Matrix<double>& pos_old, std::vector<PixelRange>& search_range_list, std::vector<Matrix<double>>& aug_img_list, double delta);

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
    Shake(
        std::vector<Matrix<int>>& orig_img_list,
        std::vector<T>& object_info, // also used as output
        double shake_width, // unit: mm
        std::vector<Camera>& cam,
        OTF& otf
    ) : _orig_img_list(orig_img_list), _object_info(object_info), _shake_width(shake_width), _cam(cam), _otf(otf) {
        // default for max_intensity
        SetMaxIntensity(255);
    };

    ~Shake() {};

    void RunShake();

    // For output: Matrix<int> residue img 
    std::vector<Matrix<int>> GetResImg() 
    {
        std::vector<Matrix<int>> res_img_list_out;

        int n_cam = _cam.size();
        for (int i = 0; i < n_cam; i ++) 
        {
            res_img_list_out.push_back(_res_img_list[i].TypeToInt());
        }

        return res_img_list_out;
    };

    void SetMaxIntensity(int max_intensity)
    {
        int n_cam = _cam.size();
        _max_intensity = std::vector<int>(n_cam, max_intensity);
    };

    void SetMaxIntensity(std::vector<int> max_intensity)
    {
        if (max_intensity.size() != _cam.size())
        {
            std::cerr << "The size of input maximum intensity vector is "
                << max_intensity.size() << "; "
                << "it does not match the size of camera" << _cam.size()
                << std::endl;
        }
        _max_intensity = max_intensity;
    };

    void SetShakeTimes (int n_loop) {
        _n_loop = n_loop;
    }

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
    }



};

#include "Shake.hpp"

#endif // !SHAKE_H

