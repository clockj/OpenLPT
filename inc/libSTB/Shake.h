#ifndef SHAKE_H
#define SHAKE_H

#include <algorithm>
#include <vector>
#include <string>
#include <omp.h>
#include <time.h>
#include <random>

#include "Matrix.h"
#include "ObjectInfo.h"
#include "Camera.h"
#include "OTF.h"


class Shake
{
public:
    // OUTPUTS //
    std::vector<Image> _imgRes_list; // residue image, size=_cam_list.n_cam_use: Original image - Particle projection image
    std::vector<double> _score_list; // score of each object: for tracer score=intensity; for bubble score=cross-correlation
    std::vector<int> _objID_remove;
    std::vector<int> _objID_keep; // _objID_keep[0] -> tr3d_list[0] object id that are kept after ghost removal
    int _n_ghost = 0;

    Shake(
        CamList const& cam_list, // all img list
        double shake_width, // unit: mm
        double score_min = 0.1, // Ghost threshold
        int n_loop = 4, // Number of shake times
        int n_thread = 0 // Number of threads
    ) : _cam_list(cam_list), _n_cam(cam_list.cam_list.size()), _n_cam_use(cam_list.useid_list.size()), _shake_width(shake_width), _score_min(score_min), _n_loop(n_loop), _n_thread(n_thread) {};

    ~Shake() {};

    // Run shake
    // if tri_only=true, only calculate residue images
    void runShake(std::vector<Tracer3D>& tr3d_list, OTF const& otf, std::vector<Image> const& imgOrig_list, bool tri_only=false);

    
private:
    // INPUTS //
    CamList const& _cam_list; 
    int _n_cam; // number of all cameras
    int _n_cam_use; // number of used cameras
    double _shake_width; // unit: mm
    double _score_min; // Ghost threshold
    int _n_loop;       // Number of shake times
    int _n_thread; // Number of threads


    //                //
    // MAIN FUNCTIONS //
    //                //

    void shakeTracers(std::vector<Tracer3D>& tr3d_list, OTF const& otf, std::vector<Image> const& imgOrig_list, bool tri_only=false);

    // Procedure for each shake
    double shakeOneTracer(Tracer3D& tr3d, OTF const& otf, double delta, double score_old);
    double shakeOneTracerGrad(Tracer3D& tr3d, OTF const& otf, double delta, double lr=1e-4);

    // Remove all tracked particles from image to get residual image.
    void calResImg(std::vector<Tracer3D> const& tr3d_list, OTF const& otf, std::vector<Image> const& imgOrig_list);

    // Remove negative pxiel and set them as zeros, this function is used to prepare residual image for the next run of IPR.
    void absResImg ();

    // Remove ghost particles.
    void removeGhost(std::vector<Tracer3D>& tr3d_list);
    void removeGhostResidue(std::vector<Tracer3D>& tr3d_list);


    //                     //
    // AUXILIARY FUNCTIONS //
    //                     //

    // Calculate the windows size for shaking
    // id: cam used id, not real cam id
    PixelRange findRegion (int id, int row, int col, int half_width_px); // a square region
    
    // Gaussian intensity for particles
    double gaussIntensity (int x, int y, Pt2D const& pt2d, std::vector<double> const& otf_param); // (x,y)=(col,row)

    // Calculate residue for shaking
    double calPointResidue (Pt3D const& pt3d, std::vector<PixelRange> const& region_list, std::vector<Image> const& imgAug_list, OTF const& otf);

    // Shaking and refine 3D position and search range
    // return final residue
    double updateTracer (Tracer3D& tr3d, std::vector<Image>& imgAug_list, std::vector<PixelRange>& region_list, OTF const& otf, double delta);
    double updateTracerGrad (Tracer3D& tr3d, std::vector<Image>& imgAug_list, std::vector<PixelRange>& region_list, OTF const& otf, double delta, double lr);

    // Update imgAug_list and region_list
    void updateImgAugList (std::vector<Image>& imgAug_list, std::vector<PixelRange>& region_list, Tracer3D const& tr3d);

    // Calculate intensity for shaken particles
    double calTracerScore (Tracer3D const& tr3d, std::vector<PixelRange> const& region_list, std::vector<Image> const& imgAug_list, OTF const& otf, double score);

};

#endif // !SHAKE_H

