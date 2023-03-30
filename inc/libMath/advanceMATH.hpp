#ifndef ADVANCEMATH_H
#define ADVANCEMATH_H

#include <vector>
#include <ctime>

#include "Matrix.h"
#include "STBCommons.hpp"
#include "Camera.h"

namespace advanceMATH
{
void Triangulation(std::vector<Camera>& cam_list, std::vector<Matrix<double>>& pt_2d_list, Matrix<double>& pt_world, double& error)
{
    Matrix<double> mtx(3, 3, 0);
    Matrix<double> pt_3d(3, 1, 0);

    std::vector<Matrix<double>> line_of_sight_list;

    for (int i = 0; i < cam_list.size(); i++)
    {
        Matrix<double> pt_2d = pt_2d_list[i];
        pt_2d = cam_list[i].ImgMMToWorld(pt_2d);
        Matrix<double> line_of_sight
            = cam_list[i].GetCenterInWorld() - pt_2d;
        line_of_sight /= line_of_sight.Mag();
        line_of_sight_list.push_back(line_of_sight);

        double nx = line_of_sight[0];
        double ny = line_of_sight[1];
        double nz = line_of_sight[2];

        Matrix<double> tmp(
            1 - nx * nx, -nx * ny, -nx * nz,
            -nx * ny, 1 - ny * ny, -ny * nz,
            -nz * nx, -nz * ny, 1 - nz * nz
        );

        pt_3d += tmp * cam_list[i].GetCenterInWorld();
        mtx += tmp;
    }

    pt_world = mtx.Inverse() * pt_3d;
    // pt_world = mtx.GaussInverse() * pt_3d;
    // pt_world = mtx.DetInverse() * pt_3d;

    // calculate errors 
    error = 0;
    for (int i = 0; i < cam_list.size(); i++)
    {
        Matrix<double> diff_vec
            = pt_world - cam_list[i].GetCenterInWorld();
        Matrix<double> line_of_sight
            = line_of_sight_list[i];
        double h = diff_vec.Mag() * diff_vec.Mag()
            - diff_vec.Dot(line_of_sight)
            * diff_vec.Dot(line_of_sight);
            
        if (h < 0 && h >= - MAGSMALLNUMBER)
        {
            h = std::fabs(h);
        }
        else if (h < - MAGSMALLNUMBER)
        {
            std::cout << "StereoMatch error: "
                << "the distance in triangulation is: "
                << "h = "
                << h
                << std::endl;
            throw;
        }

        h = std::sqrt(h);
        error += h;
    }
    error /= cam_list.size();
};


}
#endif