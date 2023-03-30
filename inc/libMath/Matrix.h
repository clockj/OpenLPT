//
//  Matrix.h
//  
//  header file of class Matrix 
//  
//  Created  by Shijie Zhong on 07/02/2022
//


#ifndef MATRIX_H
#define MATRIX_H

#include <stdlib.h>
#include <vector>
#include <iomanip>
#include <math.h>
#include <memory.h>

#include <iostream>
#include <string>
#include <sstream>
#include <fstream>

#include "myMATH.h"
#include "STBCommons.h"

template <class T>
class Matrix
{
    int _dim_x; int _dim_y;
    int _is_space = 0; // 0 for nothing installed; 1 for already claim space
    T** _mtx;
    
public:
    // Constructor
    Matrix () : Matrix(1,1) {};
    Matrix (int dim_x, int dim_y);
    Matrix (int dim_x, int dim_y, T val);
    Matrix (const Matrix<T>& mtx);
    Matrix (const Matrix<T>& vec_1, const Matrix<T>& vec_2); // generate a unit vector: v = (v_1-v_2)/|v_1-v_2|
                                                             // not applicable for int 
    Matrix (std::vector<std::vector<T>> const& mtx); // must make sure each row has the same number of columns
    Matrix (T m11, T m12,
            T m21, T m22);
    Matrix (T m11, T m12, T m13,
            T m21, T m22, T m23,
            T m31, T m32, T m33);
    Matrix (int dim_x, int dim_y, T v1, T v2, T v3); // for creating 3*1 vector
    Matrix (int dim); // generate identity square matrix
    Matrix (std::string file_name); // .csv file, _dim_x=#row, _dim_y=#col at row 1, filled with 0 if less 
    // Destructor
    ~Matrix();

    // Create/Clear space 
    void ClearMtx  ();
    void ReCreateCol (int dim_y);
    void CreateMtx (int dim_x, int dim_y, T val);
    void CreateMtx (T val);

    // Type transform
    Matrix<double> TypeToDouble ();
    Matrix<int>    TypeToInt ();

    // Get/Assign value
    T  operator() (int i, int j) const;
    T& operator() (int i, int j);
    std::vector<T> GetRow (int i);
    std::vector<T> GetCol (int j);
    double Mag ();
    double Dot (Matrix<T> const& mtx);
    // only for vector to use (return _mtx[i][0])  
    T  operator[] (int vec_i) const; 
    T& operator[] (int vec_i); 

    // sqrt(sum( xi^2 )) for all i 
    double FrobeniusNorm ();

    // Get matrix info
    int  GetDimX () const;
    int  GetDimY () const;
    void Print   ();

    // Matrix output 
    void WriteMatrix (std::string file_name);

    // Matrix calculation
    Matrix<T>& operator=  (Matrix<T> const& mtx);
    Matrix<T>  operator+  (Matrix<T> const& mtx);
    Matrix<T>& operator+= (Matrix<T> const& mtx);
    Matrix<T>  operator-  (Matrix<T> const& mtx);
    Matrix<T>& operator-= (Matrix<T> const& mtx);
    Matrix<T>  operator*  (Matrix<T> const& mtx);
    Matrix<T>& operator*= (Matrix<T> const& mtx);
    Matrix<T>  operator*  (T ratio);
    Matrix<T>& operator*= (T ratio);
    Matrix<T>  operator/  (T ratio);
    Matrix<T>& operator/= (T ratio);

    void PiecewiseMultiply (Matrix const& mtx);
    // Matrix<T> PiecewiseMultiply (Matrix& mtx);

    // Matrix manipulation
    Matrix<T> Transpose ();
    Matrix<T> Diagonal  ();
    double    Trace     ();
    Matrix<double> Inverse ();

    Matrix<double> GaussInverse ();
    void GaussForward (Matrix<double>& u_mtx, std::vector<int>& sorted_row_index);
    Matrix<double> GaussBackward (Matrix<double> b_mtx, Matrix<double> u_mtx, std::vector<int> sorted_row_index);

    Matrix<double> DetInverse  ();
};

#include "Matrix.hpp"

#endif
