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
#include <limits>
#include <algorithm>

#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <initializer_list>

#include "STBCommons.h"

template <class T>
class Matrix
{
    int _dim_row = 0; int _dim_col = 0;
    int _n = 0; // tot number of elem
    int _is_space = 0; // 0 for nothing installed; 1 for already claim space
    int mapID (int id_x, int id_y) const; // id_x*_dim_col + id_y
    T* _mtx;
    // std::vector<T> _mtx;
    
    // Create/Clear space 
    void clear  ();
    // void ReCreateCol (int dim_col);
    void create (int dim_row, int dim_col);

public:
    // Constructor
    Matrix () {};
    Matrix (Matrix<T> const& mtx); // deep copy
    Matrix (int dim_row, int dim_col, T val);

    // dim_row, dim_col must be compatible with mtx
    Matrix (std::initializer_list<std::initializer_list<T>> mtx); 

    // Load matrix from .csv file
    explicit Matrix (std::string file_name); 

    // Load matrix from input stream
    explicit Matrix (int dim_row, int dim_col, std::istream& is);
    
    // Destructor
    ~Matrix();

    // Type transform
    Matrix<double> typeToDouble ();
    Matrix<int>    typeToInt ();

    // Get/Assign value
    T  operator() (int i, int j) const;
    T& operator() (int i, int j);
    // Return _mtx[i], i = id_x*_dim_col + id_y
    // 0,1,2
    // 3,4,5 ...
    T  operator[] (int vec_i) const; 
    T& operator[] (int vec_i); 

    // Get row 
    std::vector<T> getRow(int row_i) const;
    // Get col 
    std::vector<T> getCol(int col_j) const;

    // Get matrix info
    int  getDimRow () const;
    int  getDimCol () const;
    void print (int precision = 3) const;
    const T* data() const;
    void setData (const T* data, int size);

    // Matrix output 
    void write (std::string file_name);
    void write (std::ostream& os);

    // Scalar operations
    double norm (); // sqrt(sum( xi^2 )) for all i

    // Matrix calculation
    Matrix<T>& operator=  (Matrix<T> const& mtx);
    bool       operator== (Matrix<T> const& mtx);
    bool       operator!= (Matrix<T> const& mtx);
    Matrix<T>  operator+  (Matrix<T> const& mtx) const;
    Matrix<T>& operator+= (Matrix<T> const& mtx);
    Matrix<T>  operator-  (Matrix<T> const& mtx) const;
    Matrix<T>& operator-= (Matrix<T> const& mtx);
    Matrix<T>  operator*  (Matrix<T> const& mtx) const;
    Matrix<T>& operator*= (Matrix<T> const& mtx);
    Matrix<T>  operator+  (T delta);
    Matrix<T>& operator+= (T delta);
    Matrix<T>  operator-  (T delta);
    Matrix<T>& operator-= (T delta);
    Matrix<T>  operator*  (T ratio);
    Matrix<T>& operator*= (T ratio);
    Matrix<T>  operator/  (T ratio);
    Matrix<T>& operator/= (T ratio);

    // Matrix manipulation
    Matrix<T> transpose ();

    
};

class Pt3D : public Matrix<double>
{
public:
    Pt3D () : Matrix<double>(3,1,0) {};
    Pt3D (double x, double y, double z) : Matrix<double>({{x},{y},{z}}) {};
    Pt3D (const Pt3D& pt) : Matrix<double>(pt) {};
    Pt3D (const Matrix<double>& mtx) : Matrix<double>({{mtx[0]},{mtx[1]},{mtx[2]}}) {};
    explicit Pt3D (std::string file_name) : Matrix<double>(file_name) {};
    explicit Pt3D (std::istream& is) : Matrix<double>(3,1,is) {};
};

class Pt2D : public Matrix<double>
{
public:
    Pt2D () : Matrix<double>(2,1,0) {};
    Pt2D (double x, double y) : Matrix<double>({{x},{y}}) {};
    Pt2D (const Pt2D& pt) : Matrix<double>(pt) {};
    Pt2D (const Matrix<double>& mtx) : Matrix<double>({{mtx[0]},{mtx[1]}}) {};
    explicit Pt2D (std::string file_name) : Matrix<double>(file_name) {};
    explicit Pt2D (std::istream& is) : Matrix<double>(2,1,is) {};
};

// Structure to store line
struct Line3D
{
    Pt3D pt;
    Pt3D unit_vector;
};

struct Line2D
{
    Pt2D pt;
    Pt2D unit_vector;
};

// Structure to store 3D plane
struct Plane3D
{
    Pt3D pt;
    Pt3D norm_vector;
};

// Image: matrix with double type
// Image(row_id, col_id) = intensity
// row_id = img_y, col_id = img_x
class Image : public Matrix<double>
{
public:
    Image () : Matrix<double>(1,1,0) {};
    Image (int dim_row, int dim_col, double val) : Matrix<double>(dim_row, dim_col, val) {};
    Image (std::initializer_list<std::initializer_list<double>> mtx) : Matrix<double>(mtx) {};
    Image (const Image& mtx) : Matrix<double>(mtx) {};
    Image (const Matrix<double>& mtx) : Matrix<double>(mtx) {};
    explicit Image (std::string file_name) : Matrix<double>(file_name) {};
};

#include "Matrix.hpp"

#endif
