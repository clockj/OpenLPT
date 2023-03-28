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
#include <iostream>
#include <vector>
#include <iomanip>
#include <math.h>
#include <memory.h>

#include "myMATH.hpp"

template <class T>
class Matrix
{
protected:
    int _dim_x; int _dim_y;
    int _is_space = 0; // 0 for nothing installed; 1 for already claim space
    T** _mtx;
    
public:
    // Constructor
    Matrix () : Matrix(1,1) {};
    Matrix (int dim_x, int dim_y);
    Matrix (int dim_x, int dim_y, T val);
    Matrix (const Matrix<T>& mtx);
    Matrix (std::vector<std::vector<T>> mtx);
    Matrix (T m11, T m12,
            T m21, T m22);
    Matrix (T m11, T m12, T m13,
            T m21, T m22, T m23,
            T m31, T m32, T m33);
    Matrix (int dim_x, int dim_y, T v1, T v2, T v3); // for creating 3*1 vector
    Matrix (int dim); // generate identity square matrix
    // Destructor
    ~Matrix();

    // Create/Clear space 
    void ClearMtx  ();
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
    double Dot (Matrix<T> mtx);
    // only for vector to use (return _mtx[i][0])  
    T  operator[] (int vec_i) const; 
    T& operator[] (int vec_i); 

    // sqrt(sum( xi^2 )) for all i 
    double FrobeniusNorm ();

    // Get matrix info
    std::vector<int> GetSize ();
    int  GetDimX ();
    int  GetDimY ();
    void Print   ();

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

    Matrix PiecewiseMultiply (Matrix const& mtx);
    // Matrix PiecewiseMultiply (Matrix& mtx);

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
