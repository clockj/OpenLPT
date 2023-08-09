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
	std::vector<Matrix<double>>& _orig_img_list;
	std::vector<int>& _max_intensity_list;
	std::vector<int>& _min_intensity_list;
	std::vector<Camera>& _cam_list;
	int _n_cam;
	OTF& _otf;

	bool _TRIONLY = false;
	int _n_loop_ipr = 4;
	int _n_loop_shake = 6;
	int _n_loop_ipr_reduced = 2;
	double _tol_2d = 1e-2;      // [mm], 0.25 * (x_max-x_min)/1000
	double _tol_3d = 2.4e-2;    // [mm], 0.60 * (x_max-x_min)/1000
	double _shake_width = 1e-2; // [mm]

	// OUTPUT
	std::vector<T> _object_info;
	std::vector<Matrix<double>> _res_img_list;

public:
	IPR (std::vector<Matrix<double>>& orig_img_list, std::vector<int>& max_intensity_list, std::vector<int>& min_intensity_list, std::vector<Camera>& cam_list, OTF& otf);

	~IPR () {};

	void RunIPR(std::vector<T>& object_info, bool is_reduced = false, int n_reduced = 1);
	void ReducedCamLoop (std::vector<int> cam_id, int n_cam);
	void CreateCamID (std::deque<std::vector<int>>& cam_id_all, std::vector<int> cam_id, int id, int n_rest);

	void SetTRIONLY (bool TRIONLY) {_TRIONLY = TRIONLY;};

	void SetTol2D (double tol_2d) {_tol_2d = tol_2d;};
	void SetTol3D (double tol_3d) {_tol_3d = tol_3d;};

	void SetShakeWidth (double shake_width) {_shake_width = shake_width;};
	
	void SetIPRTimes (int n_loop_ipr) { _n_loop_ipr = n_loop_ipr; };
	void SetShakeTimes (int n_loop_shake) { _n_loop_shake = n_loop_shake; };

	std::vector<Matrix<double>> GetResImgList () {return _res_img_list;};

	std::vector<T> GetObjList () {return _object_info;};
	std::vector<Matrix<double>> GetPtList () 
	{
		int n = _object_info.size();
		std::vector<Matrix<double>> pt_list(n, Matrix<double>(3,1));
		Matrix<double> pt(3,1);

		for (int i = 0; i < n; i ++)
		{
			pt_list[i] = _object_info[i].GetCenterPos();
		}

		return pt_list;
	};
};

#include "IPR.hpp"

#endif // !IPR_H
