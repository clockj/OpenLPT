//
//  Matrix.h
//  
//  header file of class Matrix 
//  
//  Created  by Shijie Zhong on 07/02/2022
//


#ifndef MATRIXBOOL_H
#define MATRIXBOOL_H

#include <stdlib.h>
#include <vector>
#include <iomanip>
#include <math.h>
#include <memory.h>

#include <iostream>
#include <string>
#include <sstream>
#include <fstream>

#include "STBCommons.h"

class Matrix_bool
{
    int _dim_x; int _dim_y;
    int _is_space = 0; // 0 for nothing installed; 1 for already claim space
    bool** _mtx;
    
public:
    // Constructor
    Matrix_bool (int dim_x, int dim_y);
    Matrix_bool (std::string file_name); // .csv file, _dim_x=#row, _dim_y=#col at row 1, filled with 0 if less 
    // Destructor
    ~Matrix_bool ();

    // Create/Clear space 
    void ClearMtx  ();
    void CreateMtx (int dim_x, int dim_y);

    // Get/Assign value
    bool  operator() (int i, int j) const;
    bool& operator() (int i, int j);
    std::vector<bool> GetRow (int i);
    std::vector<bool> GetCol (int j);
    // only for vector to use (return _mtx[i][0])  
    bool  operator[] (int vec_i) const; 
    bool& operator[] (int vec_i); 

    // Get matrix info
    int  GetDimX () const;
    int  GetDimY () const;
    void Print   ();

    // Matrix output 
    void WriteMatrix (std::string file_name);

};

#include "Matrix_bool.hpp"

#endif
