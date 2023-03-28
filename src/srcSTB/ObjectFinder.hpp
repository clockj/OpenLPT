#ifndef OBJECTFINDER_HPP
#define OBJECTFINDER_HPP

#include "ObjectFinder.h"


// Judge whether a given point is local maximum intensity or not
//  input: point on image in camera coordinate in pixel unit
//  output: bool 
template<class T>
bool ObjectFinder<T>::IsLocalMax (Matrix<int>& mtx, int i, int j)
{
    int val = mtx(i, j);
    if (
        mtx(i-1, j) > val ||
        mtx(i, j-1) > val ||
        mtx(i+1, j) > val ||
        mtx(i, j+1) > val
    ) { return false; }
    return true;
}

template<class T>
std::vector<TracerInfo> ObjectFinder<T>::FindTracer
(Matrix<int>& img, int max_intensity, int min_intensity)
{
    std::vector<TracerInfo> particle_list;

    // skip first and last rows and columns   
    for (int iter_x = 1; iter_x < img.GetDimX()-1; iter_x ++)
    {
        for (int iter_y = 1; iter_y < img.GetDimY()-1; iter_y ++)
        {
            if (img(iter_x, iter_y) > max_intensity)
            {
                // the intensity is out of range 
                std::cerr << "The intensity of pixel [" 
                          << iter_x << "," << iter_y << "] is "
                          << img(iter_x, iter_y) << ". " 
                          << "It is out of range with " 
                          << max_intensity << std::endl;
                throw;
            }

            if (img(iter_x, iter_y) >= min_intensity && IsLocalMax(img, iter_x, iter_y))
            {
                // use Gaussian distribution to find the maximum intensity
                //  => particle center position

                // Note: image row index (iter_x) => y 
                //       image col index (iter_y) => x 
                double x1 = (iter_y - 1); // +0.5;
                double x2 =  iter_y;      // +0.5;
                double x3 = (iter_y + 1); // +0.5;
                double y1 = (iter_x - 1); // +0.5;
                double y2 =  iter_x;      // +0.5;
                double y3 = (iter_x + 1); // +0.5;

                double ln_z1 = 0.0;
                double ln_z2 = 0.0;
                double ln_z3 = 0.0;

                // find the col value (coordinate: x)
                if (max_intensity == 255)
                {
                    if (img(iter_x, iter_y-1) == 0) {ln_z1 = std::log(0.0001);}
                    else {ln_z1 = myMATH::log8bit[img(iter_x, iter_y-1)];}

                    if (img(iter_x, iter_y) == 0)   {ln_z2 = std::log(0.0001);}
                    else {ln_z2 = myMATH::log8bit[img(iter_x, iter_y)];}

                    if (img(iter_x, iter_y+1) == 0) {ln_z3 = std::log(0.0001);}
                    else {ln_z3 = myMATH::log8bit[img(iter_x, iter_y+1)];}
                }
                else if (max_intensity == 65535)
                {
                    if (img(iter_x, iter_y-1) == 0) {ln_z1 = std::log(0.0001);}
                    else {ln_z1 = myMATH::log16bit[img(iter_x, iter_y-1)];}
                    if (img(iter_x, iter_y) == 0)   {ln_z2 = std::log(0.0001);}
                    else {ln_z2 = myMATH::log16bit[img(iter_x, iter_y)];}
                    if (img(iter_x, iter_y+1) == 0) {ln_z3 = std::log(0.0001);}
                    else {ln_z3 = myMATH::log16bit[img(iter_x, iter_y+1)];}
                }
                else 
                {
                    if (img(iter_x, iter_y-1) == 0) {ln_z1 = std::log(0.0001);}
                    else {ln_z1 = std::log(static_cast<double>(img(iter_x, iter_y-1)));}
                    if (img(iter_x, iter_y) == 0)   {ln_z2 = std::log(0.0001);}
                    else {ln_z2 = std::log(static_cast<double>(img(iter_x, iter_y)));}
                    if (img(iter_x, iter_y+1) == 0) {ln_z3 = std::log(0.0001);}
                    else {ln_z3 = std::log(static_cast<double>(img(iter_x,  iter_y+1)));}
                }
                double xc;
                //LONGTODO：Consider arbitrary particle size in the future.
                xc = -0.5 * (  (ln_z1 * ((x2 * x2) - (x3 * x3))) 
                             - (ln_z2 * ((x1 * x1) - (x3 * x3))) 
                             + (ln_z3 * ((x1 * x1) - (x2 * x2))) ) 
                          / (  (ln_z1 * (x3 - x2)) 
                             - (ln_z3 * (x1 - x2)) 
                             + (ln_z2 * (x1 - x3)) );
                // were these numbers valid?
                if (!std::isfinite(xc)) 
                {
                    // no -- we had a problem.  drop this particle
                    continue;
                }


                // find the row value (coordinate: y)
                if (max_intensity == 255)
                {
                    if (img(iter_x-1, iter_y) == 0) {ln_z1 = std::log(0.0001);}
                    else {ln_z1 = myMATH::log8bit[img(iter_x-1, iter_y)];}
                    if (img(iter_x+1, iter_y) == 0) {ln_z3 = std::log(0.0001);}
                    else {ln_z3 = myMATH::log8bit[img(iter_x+1, iter_y)];}
                }
                else if (max_intensity == 65535)
                {
                    if (img(iter_x-1, iter_y) == 0) {ln_z1 = std::log(0.0001);}
                    else {ln_z1 = myMATH::log16bit[img(iter_x-1, iter_y)];}

                    if (img(iter_x+1, iter_y) == 0) {ln_z3 = std::log(0.0001);}
                    else {ln_z3 = myMATH::log16bit[img(iter_x+1, iter_y)];}
                }
                else 
                {
                    if (img(iter_x-1, iter_y) == 0) {ln_z1 = std::log(0.0001);}
                    else {ln_z1 = std::log(static_cast<double>(img(iter_x-1, iter_y)));}
                    if (img(iter_x+1, iter_y) == 0) {ln_z3 = std::log(0.0001);}
                    else {ln_z3 = std::log(static_cast<double>(img(iter_x+1,  iter_y)));}
                }
                double yc; 
                yc = -0.5 * (  (ln_z1 * ((y2 * y2) - (y3 * y3))) 
                             - (ln_z2 * ((y1 * y1) - (y3 * y3))) 
                             + (ln_z3 * ((y1 * y1) - (y2 * y2))) ) 
                          / (  (ln_z1 * (y3 - y2)) 
                             - (ln_z3 * (y1 - y2)) 
                             + (ln_z2 * (y1 - y3)) );
                // were these numbers valid?
                if (!std::isfinite(yc)) 
                {
                    // no -- we had a problem.  drop this particle
                    continue;
                }


                particle_list.push_back(
                    TracerInfo(
                        Matrix<double>(
                            3,1, 
                            xc, yc, 0.0
                        )
                    )
                );

            }
        }
    }

    return particle_list;
}

template<class T>
std::vector<T> ObjectFinder<T>::FindObject(Matrix<int>& img, int max_intensity, int min_intensity)
{

    if (typeid(T) == typeid(TracerInfo))
    {
        return FindTracer
        (img, max_intensity, min_intensity);
    }
    else
    {
        std::cerr << "ObjectFinder: " 
                  << "class " << typeid(T).name()
                  << "is not included in ObjectFinder!" << std::endl;
        throw;
    }
}

#endif 