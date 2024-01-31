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

#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <initializer_list>

#include "STBCommons.h"

template <class T>
class Matrix
{
    int _dim_x; int _dim_y;
    int _n; // tot number of elem
    int _is_space = 0; // 0 for nothing installed; 1 for already claim space
    int mapID (int id_x, int id_y) const; // id_x*_dim_y + id_y
    T* _mtx;
    
    // Create/Clear space 
    void clear  ();
    // void ReCreateCol (int dim_y);
    void create (int dim_x, int dim_y);

public:
    // Constructor
    Matrix () : Matrix(1,1,0) {};
    Matrix (int dim_x, int dim_y, T val);
    // dim_x, dim_y must be compatible with mtx
    Matrix (int dim_x, int dim_y, std::initializer_list<std::initializer_list<T>> mtx); 
    Matrix (const Matrix<T>& mtx);
    // Load matrix from .csv file
    explicit Matrix (std::string file_name); 

    // Matrix (const Matrix<T>& vec_1, const Matrix<T>& vec_2); // generate a unit vector: v = (v_1-v_2)/|v_1-v_2|
    //                                                          // not applicable for int 
    
    // Destructor
    ~Matrix();

    // Type transform
    Matrix<double> typeToDouble ();
    Matrix<int>    typeToInt ();

    // Get/Assign value
    T  operator() (int i, int j) const;
    T& operator() (int i, int j);
    // Return _mtx[i], i = id_x*_dim_y + id_y
    // 0,1,2
    // 3,4,5 ...
    T  operator[] (int vec_i) const; 
    T& operator[] (int vec_i); 

    // Get row 
    std::vector<T> getRow(int row_i) const;
    // Get col 
    std::vector<T> getCol(int col_j) const;

    // Get matrix info
    int  getDimX () const;
    int  getDimY () const;
    void print   ();

    // Matrix output 
    void write (std::string file_name);

    // Scalar operations
    double norm (); // sqrt(sum( xi^2 )) for all i

    // Matrix calculation
    Matrix<T>& operator=  (Matrix<T> const& mtx);
    bool       operator== (Matrix<T> const& mtx);
    bool       operator!= (Matrix<T> const& mtx);
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

    // Matrix manipulation
    Matrix<T> transpose ();

    // // Matrix manipulation
    // Matrix<double> inverse ();

    // Matrix<double> gaussInverse ();
    // void gaussForward (Matrix<double>& u_mtx, std::vector<int>& sorted_row_index);
    // Matrix<double> gaussBackward (Matrix<double> b_mtx, Matrix<double> u_mtx, std::vector<int> sorted_row_index);

    // Matrix<double> detInverse  ();
};

class Pt3D : public Matrix<double>
{
public:
    Pt3D () : Matrix<double>(3,1,0) {};
    Pt3D (double x, double y, double z) : Matrix<double>(3,1,{{x},{y},{z}}) {};
    Pt3D (const Pt3D& pt) : Matrix<double>(pt) {};
    Pt3D (const Matrix<double>& mtx) : Matrix<double>(mtx) {};
    explicit Pt3D (std::string file_name) : Matrix<double>(file_name) {};
};

class Pt2D : public Matrix<double>
{
public:
    Pt2D () : Matrix<double>(2,1,0) {};
    Pt2D (double x, double y) : Matrix<double>(2,1,{{x},{y}}) {};
    Pt2D (const Pt2D& pt) : Matrix<double>(pt) {};
    Pt2D (const Matrix<double>& mtx) : Matrix<double>(mtx) {};
    explicit Pt2D (std::string file_name) : Matrix<double>(file_name) {};
};

// Structure to store line
struct Line3D
{
    Pt3D _pt;
    Pt3D _unit_vector;
};

struct Line2D
{
    Pt2D _pt;
    Pt2D _unit_vector;
};

class Image : public Matrix<double>
{
public:
    Image () : Matrix<double>(1,1,0) {};
    Image (int dim_x, int dim_y, double val) : Matrix<double>(dim_x, dim_y, val) {};
    Image (int dim_x, int dim_y, std::initializer_list<std::initializer_list<double>> mtx) : Matrix<double>(dim_x, dim_y, mtx) {};
    Image (const Image& mtx) : Matrix<double>(mtx) {};
    Image (const Matrix<double>& mtx) : Matrix<double>(mtx) {};
    explicit Image (std::string file_name) : Matrix<double>(file_name) {};
};

#include "Matrix.hpp"

#endif
