#ifndef KALMAN_H
#define KALMAN_H

#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>

#include "STBCommons.h"
#include "myMATH.h"
#include "Matrix.h"

class KalmanFilter
{
public:
    Matrix<double> _F;  // State transition model
    Matrix<double> _H;  // Measurement model
    Matrix<double> _Q;  // Process noise covariance
    Matrix<double> _R;  // Measurement noise covariance
    Matrix<double> _x;  // Initial state estimate
    Matrix<double> _P;  // Initial covariance estimate

    // Functions //
    KalmanFilter() {};
    KalmanFilter(const Matrix<double>& F, const Matrix<double>& H, const Matrix<double>& Q, const Matrix<double>& R, const Matrix<double>& x0, const Matrix<double>& P0);
    KalmanFilter(KalmanFilter const& kf);
    ~KalmanFilter() {};

    // Predict the next position of the track
    // but not update the obj2d list
    void predict();

    // Update the filter with the new observation
    void update(Matrix<double> z);

};

#endif