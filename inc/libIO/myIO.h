#ifndef MYIO_H
#define MYIO_H

#include <iostream>
#include <string>
#include <sstream>
#include <fstream>

#include "Matrix.h"
#include "ObjectInfo.h"

namespace myIO
{

template<class T>
void WriteMatrix (std::string file_name, std::vector<std::vector<T>> const& mtx);

void WriteTracerPos (std::string file_name, std::vector<TracerInfo> const& object_list);

// .csv file 
// not recommended for large matrix
template<class T>
void LoadMatrix (std::string file_name, std::vector<std::vector<T>>& mtx);

}


#endif