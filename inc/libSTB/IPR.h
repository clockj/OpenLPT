#ifndef IPR_H
#define IPR_H

#include <vector>
#include "Matrix.h"
#include "Camera.h"
#include "StereoMatch.h"
#include "Shake.h"
#include "OTF.h"
#include "ObjectFinder.h"

template<class T>
class IPR 
{
protected:
	// INPUT
	std::vector<Matrix<int>>& _orig_img_list;
	std::vector<int>& _max_intensity_list;
	std::vector<int>& _min_intensity_list;
	std::vector<Camera>& _cam_list;
	OTF& _otf;

	int _n_loop_ipr = 2;
	int _n_loop_shake = 6;
	double _tol_2d = 2e-2; // 3e-2 [mm]
	double _tol_3d = 5e-1; // 8e-2 [mm]
	double _shake_width = 1e-1; // [mm]

	// OUTPUT
	std::vector<T> _object_info;

public:
	IPR (std::vector<Matrix<int>>& orig_img_list, std::vector<int>& max_intensity_list, std::vector<int>& min_intensity_list, std::vector<Camera>& cam_list, OTF& otf);

	~IPR () {};

	std::vector<T> RunIPR();

	void SetTol2D (double tol_2d) { _tol_2d = tol_2d; };
	void SetTol3D (double tol_3d) { _tol_3d = tol_3d; };

	void SetShakeWidth (double shake_width) {_shake_width = shake_width;};
	
	void SetIPRTimes (int n_loop_ipr) { _n_loop_ipr = n_loop_ipr; };
	void SetShakeTimes (int n_loop_shake) { _n_loop_shake = n_loop_shake; };
};

#include "IPR.hpp"

#endif // !IPR_H