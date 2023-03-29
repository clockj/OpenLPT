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

#include "STBCommons.h"

namespace myMATH
{

// Sorting 
template<class T>
void Swap (std::vector<T>& nums, int i, int j);
// return bubble sorted index (small->large)
template<class T>
void BubbleSort (std::vector<int>& sorted_id_list, std::vector<T> const& nums);

template<class T>
void Merge (std::vector<int>& B, int id_begin, int id_middle, int id_end, std::vector<int>& A, std::vector<T> const& C);
template<class T>
void SplitMerge (std::vector<int>& B, int id_begin, int id_end, std::vector<int>& A, std::vector<T> const& C);
// return merge sorted index (small->large)
template<class T>
void MergeSort (std::vector<int>& sorted_id_list, std::vector<T> const& nums);


// Max or Min 
template<class T>
T Max(std::vector<T> const& nums, SortTypeID type_id = type_ms);
template<class T>
T Min(std::vector<T> const& nums, SortTypeID type_id = type_ms);


// Calculate median
template<class T>
double Median (std::vector<T> const& nums);


// Find outlier in a vector
template<class T>
void IsOutlier (std::vector<bool>& judge, std::vector<T> const& nums);


// Linespace
// n >= 2, including min and max, n_gap = n-1
std::vector<double> Linspace (double min, double max, int n);

double TriLinearInterp(AxisLimit& grid_limit, std::vector<double>& value, std::vector<double>& pt_vec);


} 


#endif