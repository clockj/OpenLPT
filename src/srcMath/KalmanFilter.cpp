#include "KalmanFilter.h"

KalmanFilter::KalmanFilter(const KalmanFilter& kf)
    : _F(kf._F), _H(kf._H), _Q(kf._Q), _R(kf._R), _x(kf._x), _P(kf._P) {}

KalmanFilter::KalmanFilter(const Matrix<double>& F, const Matrix<double>& H, const Matrix<double>& Q, const Matrix<double>& R, const Matrix<double>& x0, const Matrix<double>& P0)
    : _F(F), _H(H), _Q(Q), _R(R), _x(x0), _P(P0) {}

void KalmanFilter::predict()
{
    // Predict the next state and update the covariance matrix.
    _x = _F * _x;
    _P = _F * _P * _F.transpose() + _Q;
}

void KalmanFilter::update(Matrix<double> z)
{
    // Update the state estimate based on the measurement z.
    z -= _H * _x;
    Matrix<double> S = _H * _P * _H.transpose() + _R;
    Matrix<double> K = _P * _H.transpose() * myMATH::inverse(S);

    _x = _x + K * z;
    _P = (myMATH::eye<double>(_P.getDimRow()) - K * _H) * _P;
}