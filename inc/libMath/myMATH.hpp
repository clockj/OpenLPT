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

#include "STBCommons.hpp"

namespace myMATH
{

// Sorting 
template<class T>
void Swap (std::vector<T>& nums, int i, int j)
{
    T temp = nums[j];
    nums[j] = nums[i];
    nums[i] = temp;
};
// return bubble sorted index (small->large)
template<class T>
void BubbleSort (std::vector<int>& sorted_id_list, std::vector<T> const& nums)
{
    if (sorted_id_list.size() != nums.size())
    {
        std::cerr << "myMATH::BubbleSort: size unequal" << std::endl;
        throw -1;
    }

    int size = nums.size();
    std::vector<T> nums_copy = nums;

    for (int i = 0; i < size; i ++)
    {
        sorted_id_list[i] = i;
    }

    for (int i = 0; i < size-1; i++)
    {
        for (int j = 0; j < size-1-i; j++)
        {
            if (nums_copy[j] > nums_copy[j+1])
            {   
                Swap<T>  (nums_copy, j, j+1);
                Swap<int>(sorted_id_list, j, j+1);
            }
        }
    }
};

template<class T>
void Merge (std::vector<int>& B, int id_begin, int id_middle, int id_end, std::vector<int>& A, std::vector<T> const& C)
{
    int i = id_begin;
    int j = id_middle;

    for (int k = id_begin; k < id_end; k ++)
    {
        if (i < id_middle && j >= id_end)
        {
            A[k] = B[i];
            i ++;
        }
        else if (i < id_middle && C[B[i]] < C[B[j]])
        {
            A[k] = B[i];
            i ++;
        }
        else
        {
            A[k] = B[j];
            j ++; 
        }
    }
};
template<class T>
void SplitMerge (std::vector<int>& B, int id_begin, int id_end, std::vector<int>& A, std::vector<T> const& C)
{
    if (id_end - id_begin < 2)
    {
        return;
    }
    int id_middle = (id_begin + id_end) / 2;
    SplitMerge(A, id_begin, id_middle, B, C);
    SplitMerge(A, id_middle, id_end, B, C);
    Merge(B, id_begin, id_middle, id_end, A, C);
};
// return merge sorted index (small->large)
template<class T>
void MergeSort (std::vector<int>& sorted_id_list, std::vector<T> const& nums)
{
    int n = nums.size();

    std::vector<int> B(n);
    for (int i = 0; i < n; i ++)
    {
        sorted_id_list[i] = i;
        B[i] = i;
    }

    SplitMerge<T> (B, 0, n, sorted_id_list, nums);
};


// Max or Min 
template<class T>
T Max(std::vector<T> const& nums, SortTypeID type_id = type_ms)
{
    int size = nums.size();
    if (size < 1)
    {
        std::cerr << "myMATH::Max: size=" << size << std::endl;
        throw -1;
    }

    std::vector<int> sorted_id_list(size);
    switch (type_id)
    {
    case type_ms:  
        MergeSort<T> (sorted_id_list, nums);
        return nums[sorted_id_list[size-1]];
    case type_bs:
        BubbleSort<T> (sorted_id_list, nums);
        return nums[sorted_id_list[size-1]];
    default:
        std::cerr << "myMATH::Max: no such sorting" << std::endl;
        throw -2;
    }
};
template<class T>
T Min(std::vector<T> const& nums, SortTypeID type_id = type_ms)
{
    int size = nums.size();
    std::vector<int> sorted_id_list(size);
    switch (type_id)
    {
    case type_ms:  
        MergeSort<T> (sorted_id_list, nums);
        return nums[sorted_id_list[0]];
    case type_bs:
        BubbleSort<T> (sorted_id_list, nums);
        return nums[sorted_id_list[0]];
    default:
        std::cerr << "myMATH::Max: no such type" << std::endl;
        return 0;
    }
};


// Calculate median
template<class T>
double Median (std::vector<T> const& nums)
{
    int n = nums.size();
    std::vector<int> sort_index(n);
    MergeSort<T> (sort_index, nums);

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
void IsOutlier (std::vector<bool>& judge, std::vector<T> const& nums)
{
    int n = nums.size();
    if (judge.size() != n)
    {
        std::cerr << "myMATH::IsOutlier: size unequal" << std::endl;
        throw -1;
    }

    double median = Median<T>(nums);
    std::vector<double> deviation(n, 0);

    for (int i = 0; i < n; i ++)
    {
        deviation[i] = std::fabs(nums[i] - median);
    }

    double amad = 1.4826022185056018 * Median(deviation);
    double lb = median - 3 * amad;
    double rb = median + 3 * amad;

    for (int i = 0; i < n; i ++)
    {
        judge[i] = (nums[i] < lb) || (nums[i] > rb);
    }
};


// Linespace
// n >= 2, including min and max, n_gap = n-1
std::vector<double> Linspace (double min, double max, int n)
{
    // both min and max are stored in the vector
    if (n < 2)
    {
        std::cerr << "myMATH::Linspace: Cannot generate line space: "
                  << "n = " << n << ", "
                  << "less than 2!"
                  << std::endl;
        throw -1;
    }

    double delta = (max - min) / (n - 1);

    std::vector<double> res(n);
    for (int i = 0; i < n; i ++)
    {
        res[n] = (min + delta * i);
    }

    return res;
};

double TriLinearInterp(AxisLimit& grid_limit, std::vector<double>& value, std::vector<double>& pt_vec)
{
    // [https://en.wikipedia.org/wiki/Trilinear_interpolation]
    // value: 8
    //        c_000 [x_min, y_min, z_min]
    //        c_100 [x_max, y_min, z_min]
    //        c_101 [x_max, y_min, z_max]
    //        c_001 [x_min, y_min, z_max]
    //        c_010 [x_min, y_max, z_min]
    //        c_110 [x_max, y_max, z_min]
    //        c_111 [x_max, y_max, z_max]
    //        c_011 [x_min, y_max, z_max]
    // pt_world: 3 
    //           x, y, z
    
    // Key behavior:
    // If pt_world is outside the limit:
    //  set it as the value on the boundary 

    double x_0 = grid_limit._x_min;
    double x_1 = grid_limit._x_max;
    double y_0 = grid_limit._y_min;
    double y_1 = grid_limit._y_max;
    double z_0 = grid_limit._z_min;
    double z_1 = grid_limit._z_max;

    double c_000 = value[0];
    double c_100 = value[1];
    double c_101 = value[2];
    double c_001 = value[3];
    double c_010 = value[4];
    double c_110 = value[5];
    double c_111 = value[6];
    double c_011 = value[7];

    double x = pt_vec[0];
    double y = pt_vec[1];
    double z = pt_vec[2];

    double x_d = (x - x_0) / (x_1 - x_0);
    double y_d = (y - y_0) / (y_1 - y_0);
    double z_d = (z - z_0) / (z_1 - z_0);

    if (x_d > 1 || x_d < 0 || 
        y_d > 1 || y_d < 0 || 
        z_d > 1 || z_d < 0)
    {
        std::cerr << "myMATH::TriLinearInterp error: out of range" 
                  << "(x,y,z) = (" 
                  << x << ","
                  << y << ","
                  << z << ")"
                  << std::endl;
        throw -3;
    }

    double c_00 = c_000 * (1 - x_d) + c_100 * x_d;
    double c_01 = c_001 * (1 - x_d) + c_101 * x_d;
    double c_10 = c_010 * (1 - x_d) + c_110 * x_d;
    double c_11 = c_011 * (1 - x_d) + c_111 * x_d;

    double c_0 = c_00 * (1 - y_d) + c_10 * y_d;
    double c_1 = c_01 * (1 - y_d) + c_11 * y_d;

    double c = c_0 * (1 - z_d) + c_1 * z_d;

    return c;
};


} 

#endif