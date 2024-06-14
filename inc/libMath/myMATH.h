//
//  MATH.h
//  
//  header file of class Matrix 
//  
//  Created  by Shijie Zhong on 07/02/2022
//


#ifndef MYMATH_H
#define MYMATH_H

#include <cstdlib>
#include <vector>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <typeinfo>
#include <string>

#include "STBCommons.h"
#include "Matrix.h"

namespace myMATH
{

// Sort index
template<class T>
void sortID (std::vector<int>& sort_index, std::vector<T> const& nums)
{
    int n = nums.size();
    if (sort_index.size() != n)
    {
        sort_index.resize(n);
    }

    // Initialize sort_index
    for (int i = 0; i < n; i ++)
    {
        sort_index[i] = i;
    }

    // Custom comparator function for sorting indices based on corresponding values
    auto comparator = [nums](size_t i, size_t j) {
        return nums[i] < nums[j];
    };

    // Use std::sort with custom comparator
    std::sort(sort_index.begin(), sort_index.end(), comparator);
};


// Calculate median
template<class T>
double getMedian (std::vector<T> const& nums)
{
    int n = nums.size();
    std::vector<int> sort_index(n);
    sortID<T> (sort_index, nums);

    if (n & 1) // odd number: check if the last digit is 1
    {
        // median is the middle number 
        return double(nums[sort_index[(n+1)/2-1]]);
    }
    else // even number 
    {
        return double((nums[sort_index[n/2-1]] + nums[sort_index[n/2]]) / 2);
    }
};


// Find outlier in a vector
template<class T>
void isOutlier (std::vector<int>& judge, std::vector<T> const& nums)
{
    int n = nums.size();
    if (judge.size() != n)
    {
        judge.resize(n);
    }

    double median = getMedian<T> (nums);
    std::vector<double> deviation(n, 0);

    for (int i = 0; i < n; i ++)
    {
        deviation[i] = std::fabs(nums[i] - median);
    }

    double amad = 1.4826022185056018 * getMedian<double> (deviation);
    double lb = median - 3 * amad;
    double rb = median + 3 * amad;

    for (int i = 0; i < n; i ++)
    {
        if ((nums[i] < lb) || (nums[i] > rb))
        {
            judge[i] = 1;
        }
        else
        {
            judge[i] = 0;
        }
    }
};


// Linespace
// n >= 2, including min and max, n_gap = n-1
std::vector<double> linspace (double min, double max, int n);


// Trilinear interpolation
double triLinearInterp(AxisLimit const& grid_limit, std::vector<double> const& value, std::vector<double> const& pt_vec);


// Create unit vector 
Pt3D createUnitVector (Pt3D const& pt1, Pt3D const& pt2);
Pt2D createUnitVector (Pt2D const& pt1, Pt2D const& pt2);


// Calculate dot product
double dot (Pt3D const& pt1, Pt3D const& pt2);
double dot (Pt2D const& pt1, Pt2D const& pt2);

// TODO: add dist2
// Calculate the distance square between two points
double dist2 (Pt3D const& pt1, Pt3D const& pt2);
double dist2 (Pt2D const& pt1, Pt2D const& pt2);


// Calculate the distance square between point and line  
double dist2 (Pt3D const& pt, Line3D const& line);
double dist2 (Pt2D const& pt, Line2D const& line);


// Calculate the distance between two points
double dist (Pt3D const& pt1, Pt3D const& pt2);
double dist (Pt2D const& pt1, Pt2D const& pt2);


// Calculate the distance between point and line  
double dist (Pt3D const& pt, Line3D const& line);
double dist (Pt2D const& pt, Line2D const& line);


// Triangulation
void triangulation (Pt3D& pt_world, double& error,
                    std::vector<Line3D> const& line_of_sight_list);


// Find cross points of two 2d lines
Pt2D crossPoint (Line2D const& line1, Line2D const& line2);


// Create identity matrix 
template<class T>
Matrix<T> eye (int n)
{
    Matrix<T> res(n,n,0);
    for (int i = 0; i < n; i ++)
    {
        res(i,i) = 1;
    }
    return res;
};


// Calculate piecewise product
template<class T>
Matrix<T> piecewiseProduct (Matrix<T> const& mtx1, Matrix<T> const& mtx2)
{
    if ((mtx1.getDimRow() != mtx2.getDimRow()) || (mtx1.getDimCol() != mtx2.getDimCol()))
    {
        std::cerr << "myMATH::PiecewiseProduct: size unequal" << std::endl;
        throw error_size;
    }

    Matrix<T> res(mtx1);
    for (int i = 0; i < mtx1.getDimRow(); i ++)
    {
        for (int j = 0; j < mtx1.getDimCol(); j ++)
        {
            res(i,j) *= mtx2(i,j);
        }
    }

    return res;
};


// Internal functions
namespace internal
{
    template<class T> 
    void gaussForward(Matrix<T>& u_mtx, std::vector<int>& sortedRowIndex)
    {
        // Solve Amat * Xmat = Bmat,
        //  where Amat (n*n), Xmat (n*m), Bmat (n*m).
        // In GaussForward, Amat = Lmat * u_mtx, 
        //  we are trying to solve for u_mtx (upper triangle matrix).
        // The actions (Lmat^-1) are stored in 
        //  the lower triangle part of u_mtx (multipliers),
        //  and the sortedRowIndex.

        int n = u_mtx.getDimRow(); // dimension of the matrix: n * n

        std::vector<T> row(n, 0.0);
        std::vector<T> maxOfrow(n, 0.0);

        for (int iter_x = 0; iter_x < n; iter_x ++)
        {
            sortedRowIndex[iter_x] = iter_x;
            row = u_mtx.getRow(iter_x);
            maxOfrow[iter_x] = *std::max_element(row.begin(), row.end());
        }

        for (int iter_y = 0; iter_y < n-1; iter_y ++)
        {
            T maxPivot = 0; // "biggest pivot" of each col
            int t = 0;      // row index of "biggest pivot"
            T pivot;        // "relative" pivot value
            for (int iter_x = iter_y; iter_x < n; iter_x ++)
            {
                pivot = std::fabs(u_mtx(sortedRowIndex[iter_x], iter_y) 
                                / maxOfrow[sortedRowIndex[iter_x]]);
                if (pivot > maxPivot)
                {
                    maxPivot = pivot;
                    t = iter_x;
                }
            }
            std::swap(sortedRowIndex[t], sortedRowIndex[iter_y]); // move the row with "biggest pivot" 
                                                                  // to the first row of the sub-matrix

            T m; // multiplier during elimination
            for (int iter_x = iter_y+1; iter_x < n; iter_x ++)
            {
                m = u_mtx(sortedRowIndex[iter_x], iter_y) 
                    / u_mtx(sortedRowIndex[iter_y], iter_y);
                u_mtx(sortedRowIndex[iter_x], iter_y) = m;  // it should be 0. 
                                                            // here is only used to store the multiplier m

                for (int iter = iter_y+1; iter < n; iter ++)
                {
                    u_mtx(sortedRowIndex[iter_x], iter) -= m * u_mtx(sortedRowIndex[iter_y], iter);
                }
            }
        }
    };

    template<class T> 
    Matrix<T> gaussBackward(Matrix<T>& b_mtx, Matrix<T>& u_mtx, std::vector<int>& sortedRowIndex)
    {
        // Solve Amat * Xmat = Bmat,
        //  where Amat (n*n), Xmat (n*m), Bmat (n*m).
        // In GaussForward, Amat = Lmat * u_mtx, 
        //  where Lmat (n*n) and u_mtx (n*n),
        //  we are trying to solve for u_mtx (upper triangle matrix).
        // The actions (Lmat^-1) are stored in 
        //  the lower triangle part of u_mtx (multipliers),
        //  and the sortedRowIndex.
        // In GaussBackward, we are trying to solve for Xmat,
        //  where Xmat = u_mtx^(-1) * Lmat^(-1) * Bmat

        int n = u_mtx.getDimRow(); // dimension of the matrix: n * n
        Matrix<T> res(b_mtx);
        T z; // right-hand-side variable 

        for (int iter = 0; iter < n; iter ++)
        {
            for (int iter_y = 0; iter_y < n-1; iter_y ++)
            {
                for (int iter_x = iter_y+1; iter_x < n; iter_x ++)
                {
                    b_mtx(sortedRowIndex[iter_x], iter) -= (u_mtx(sortedRowIndex[iter_x], iter_y) 
                                                          * b_mtx(sortedRowIndex[iter_y], iter));
                }
            }

            // sol to the last variable
            res(n-1, iter) = b_mtx(sortedRowIndex[n-1], iter) 
                           / u_mtx(sortedRowIndex[n-1], n-1); 
            
            for (int iter_x = n-2; iter_x > -1; iter_x --)
            {
                z = b_mtx(sortedRowIndex[iter_x], iter);
                for (int iter_y = iter_x+1; iter_y < n; iter_y ++)
                {
                    z -= (u_mtx(sortedRowIndex[iter_x], iter_y) * res(iter_y, iter));
                }
                res(iter_x, iter) = z / u_mtx(sortedRowIndex[iter_x], iter_x);
            }
            
        }

        return res;
    };

    template<class T> 
    Matrix<T> gaussInverse (Matrix<T> const& mtx)
    {
        // Gaussian Elimination
        int n = mtx.getDimRow();
        std::vector<int> sortedRowIndex(n);
        Matrix<T> u_mtx (mtx); // deep copy
        Matrix<T> b_mtx(eye<T>(n)); // generate identity matrix

        gaussForward(u_mtx, sortedRowIndex);
        Matrix<T> res = gaussBackward(b_mtx, u_mtx, sortedRowIndex);

        return res;
    };

    template<class T> 
    Matrix<T> detInverse (Matrix<T> const& mtx)
    {
        if (mtx.getDimRow() != mtx.getDimCol())
        {
            std::cerr << "The matrix is not a square matrix!"
                    << "There is no inverse."
                    << std::endl;
            throw error_size;
        }

        int n = mtx.getDimRow();
        if (n == 2)
        {
            Matrix<T> res(2,2,0);
            //        [0]    [1]
            // [0] res[0] res[1]
            // [1] res[2] res[3]

            T det = mtx[0] * mtx[3] - mtx[1] * mtx[2];
            if (std::fabs(det) < SMALLNUMBER)
            {
                std::cout << "myMATH::detInverse warning at line " 
                          << __LINE__ << ": "
                          << "The determinant is too small, "
                          << det << ". "
                          << "Singular matrix!"
                          << std::endl;
                det = SMALLNUMBER;
                mtx.print();
            }

            res[0] =   mtx[3] / det;
            res[1] = - mtx[1] / det;
            res[2] = - mtx[2] / det;
            res[3] =   mtx[0] / det;

            return res;
        }
        else if (n == 3)
        {
            Matrix<T> res(3,3,0);
            //        [0]    [1]    [2]
            // [0] res[0] res[1] res[2]
            // [1] res[3] res[4] res[5]
            // [2] res[6] res[7] res[8]

            T det =   mtx[0] * ( mtx[4]*mtx[8] - mtx[5]*mtx[7] )
                    - mtx[1] * ( mtx[3]*mtx[8] - mtx[5]*mtx[6] ) 
                    + mtx[2] * ( mtx[3]*mtx[7] - mtx[4]*mtx[6] );
            if (std::fabs(det) < SMALLNUMBER)
            {
                std::cout << "myMATH::detInverse warning at line " 
                          << __LINE__ << ": "
                          << "The determinant is too small, "
                          << det << ". "
                          << "Singular matrix!"
                          << std::endl;
            }
            
            res[0] =   (mtx[4]*mtx[8] - mtx[5]*mtx[7]) / det;
            res[3] = - (mtx[3]*mtx[8] - mtx[5]*mtx[6]) / det;
            res[6] =   (mtx[3]*mtx[7] - mtx[4]*mtx[6]) / det;

            res[1] = - (mtx[1]*mtx[8] - mtx[2]*mtx[7]) / det;
            res[4] =   (mtx[0]*mtx[8] - mtx[2]*mtx[6]) / det;
            res[7] = - (mtx[0]*mtx[7] - mtx[1]*mtx[6]) / det;

            res[2] =   (mtx[1]*mtx[5] - mtx[2]*mtx[4]) / det;
            res[5] = - (mtx[0]*mtx[5] - mtx[2]*mtx[3]) / det;
            res[8] =   (mtx[0]*mtx[4] - mtx[1]*mtx[3]) / det;

            return res;
        }
        else
        {
            std::cerr << "DetInverse only supports for 2x2 or 3x3 matrix!"
                    << std::endl;
            throw error_size;
        }
        
    };

};


// Calculate matrix inverse 
template<class T> 
Matrix<T> inverse (Matrix<T> const& mtx, std::string method="gauss")
{   
    std::string type = typeid(T).name();
    if (type.find("int") != std::string::npos)
    {
        std::cerr << "myMATH::Inverse: not supported for int" << std::endl;
        throw error_type;
    }
    if (mtx.getDimRow() != mtx.getDimCol())
    {
        std::cerr << "myMATH::Inverse: not a square matrix" << std::endl;
        throw error_size;
    }

    if (method == "gauss")
    {
        return internal::gaussInverse<T>(mtx);
    }
    else if (method == "det")
    {
        return internal::detInverse<T>(mtx);
    }
    else
    {
        std::cerr << "myMATH::Inverse: method not supported" << std::endl;
        throw error_type;
    }
};


// Calculate trace of a matrix
template<class T>
T trace (Matrix<T> const& mtx)
{
    int dim_row = mtx.getDimRow();
    int dim_col = mtx.getDimCol();

    if (dim_row != dim_col)
    {
        std::cerr << "myMATH::trace: matrix is not square" << std::endl;
        throw error_size;
    }

    T res = 0;   
    for (int i = 0; i < dim_row; i ++)
    {
        res += mtx(i,i);
    }
    return res;
};


// Check local maximum
template<class T>
bool isLocalMax (Matrix<T> const& mtx, int row_id, int col_id)
{
    const int n_row = mtx.getDimRow();
    const int n_col = mtx.getDimCol();
    if (row_id < 0 || row_id >= n_row || col_id < 0 || col_id >= n_col)
    {
        std::cerr << "myMATH::IsLocalMax: index out of range" << std::endl;
        std::cerr << "row_id: " << row_id << ", col_id: " << col_id << std::endl;
        std::cerr << "n_row: " << n_row << ", n_col: " << n_col << std::endl;
        throw error_range;
    }

    const int move_row[4] = {-1, 0, 1, 0};
    const int move_col[4] = {0, -1, 0, 1};

    T val = mtx(row_id, col_id);
    for (int i = 0; i < 4; i ++)
    {
        if (row_id + move_row[i] < 0 || row_id + move_row[i] >= n_row ||
            col_id + move_col[i] < 0 || col_id + move_col[i] >= n_col)
        {
            continue;
        }

        if (mtx(row_id+move_row[i], col_id+move_col[i]) > val)
        {
            return false;
        }
    }

    return true;
};


// Fit quadratic function
// y = coeff[0] + coeff[1] * x + coeff[2] * x^2 + ...
void polyfit (std::vector<double>& coeff, std::vector<double> const& x, std::vector<double> const& y, int order);

}



#endif