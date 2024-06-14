//
//	ObjectFinder.h 
//
//	Base class for finding object from images
//
//	Created by Shijie Zhong 07/06/2022 
//

#ifndef OBJECTFINDER_H
#define OBJECTFINDER_H

#include <vector>
#include <typeinfo>

#include "ObjectInfo.h"
#include "Matrix.h"
#include "STBCommons.h"
#include "myMATH.h"


class ObjectFinder2D
{
private:
    void findTracer2D(std::vector<Tracer2D>& tr2d_list, Image const& img, int max_intensity, int min_intensity, double r_px=2);

    void findTracer2D(std::vector<Tracer2D>& tr2d_list, Image const& img, int max_intensity, int min_intensity, double r_px, PixelRange const& region);

public:
    ObjectFinder2D() {};
    ~ObjectFinder2D() {};

    // Find object position
    //  input: intensity matrix, maximum intensity (2^bit_per_sample-1)
    //  output: a vector with all the particles positions （x_pixel(col_id), y_pixel(row_id))
    template <class T> 
    void findObject2D(std::vector<T>& obj2d_list, Image const& img, std::vector<double> const& properties);

    template <class T>
    void findObject2D(std::vector<T>& obj2d_list, Image const& img, std::vector<double> const& properties, PixelRange const& region);

};

#include "ObjectFinder.hpp"

#endif