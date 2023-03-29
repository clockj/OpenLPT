#ifndef MYIO_H
#define MYIO_H

#include <iostream>
#include <string>
#include <sstream>
#include <fstream>

#include "Matrix.h"

namespace myIO
{

template<class T>
void WriteMatrix (std::string file_name, std::vector<std::vector<T>> const& mtx);

// .csv file 
// not recommended for large matrix
template<class T>
void LoadMatrix (std::string file_name, std::vector<std::vector<T>>& mtx);

}


#endif