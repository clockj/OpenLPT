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

// template<ObjectType T>
template <class T> 
class ObjectFinder
{
protected:
    bool IsLocalMax (Matrix<int> const& mtx, int i, int j);
    std::vector<TracerInfo> FindTracer(Matrix<int> const& img, int max_intensity, int min_intensity);
public:
    ObjectFinder() {};
    virtual ~ObjectFinder() {};

    // Find object position
    //  input: intensity matrix, maximum intensity (2^bit_per_sample-1)
    //  output: a vector with all the particles positions
    std::vector<T> FindObject(Matrix<int> const& img, int max_intensity, int min_intensity);

};

#include "ObjectFinder.hpp"

#endif