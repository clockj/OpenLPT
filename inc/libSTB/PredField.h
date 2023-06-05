#ifndef PREDFIELD_H
#define PREDFIELD_H

#include "Matrix.h"
#include "myMATH.h"

class PredField
{
private:
    // X,Y,Z limits of view area
    AxisLimit _limit;
    // # of grids in X,Y,Z
    std::vector<int> _n_xyz;
    int _n_tot;
    // grid spacing in X,Y,Z
    std::vector<double> _grid_size;
    
    // grid points 
    // i in 0 ~ n_grid: i = index_x*_ny*_nz + index_y*_nz + index_z
    // for index_x in 0:_nx-1 
    //     for index_y in 0:_ny-1 
    //         for index_z in 0:_nz-1  
    //             i ++
    Matrix<double> _grid;

    // 3D particles at prev and curr frame
    std::vector<Matrix<double>>& _pt_list_prev;
    std::vector<Matrix<double>>& _pt_list_curr;

    // radius of interrogation sphere
    double _r;
    // parameters for displacement map
    int _disp_map_res = 10; 
    int _size;
    double _m, _c;

    Matrix<double> _disp_field;

    // Set grid
    void SetGrid();

public:
    // constructor : To get predictive field( takes the previous frame 3D pos, field parameters from file, frame, matlabFlag )
    PredField(AxisLimit& limit, std::vector<int>& n_xyz, std::vector<Matrix<double>>& pt_list_prev, std::vector<Matrix<double>>& pt_list_curr, double r) 
        : _limit(limit), _n_xyz(n_xyz), _n_tot(n_xyz[0]*n_xyz[1]*n_xyz[2]), _grid(3, _n_tot), _pt_list_prev(pt_list_prev), _pt_list_curr(pt_list_curr), _r(r)
    {
        SetGrid();

        // converting the Map index (i) to displacement (dx) as i = m*dx + c;
        _size = 4 * _disp_map_res * _r + 1; // dispMapRes = 10
        _m = (_size - 1) / (4 * _r);
        _c = (_size - 1) / 2;
    };

    // destructor
    ~PredField() {};

    int MapGridIndexToID (int index_x, int index_y, int index_z)
    {
        return index_x*_n_xyz[1]*_n_xyz[2] + index_y*_n_xyz[2] + index_z;
    }

    void Field();
    void FindVolPt(std::vector<int>& pt_list_id, double rsqr, int grid, FrameTypeID frame);
    void DispMap(std::vector<std::vector<std::vector<int>>>& dispMap, std::vector<std::vector<double>>& displacements);
    std::vector<double> DispMapPeak(std::vector<std::vector<std::vector<int>>>& dispMap);
    std::vector<int> MaxElemD(double*** &array);
    double Gauss1DPeak(double y1, double v1, double y2, double v2, double y3, double v3);

    // getting the necessary variables
    Matrix<double> GetGrid() {
        return _grid;
    }
    Matrix<double> GetField() {
        return _disp_field;
    }
    std::vector<Matrix<double>> GetCurrPos3D() {
        return _pt_list_curr;
    }
    std::vector<Matrix<double>> GetPrevPos3D() {
        return _pt_list_prev;
    }

    Matrix<double> ParticleInterpolation(Matrix<double> const& pos3D);

    // void MatfileSave(vector<double*> pos, string name);
    // void MatfileSave(vector<double> pos[3], string name);
    // void MatfileSave(vector< vector<double> > pos, string name);
    // void SaveField(string file_path);
    // void Load_field(string path);

};


#endif