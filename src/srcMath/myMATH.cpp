#include "myMATH.h"

namespace myMATH
{

// Linespace
// n >= 2, including min and max, n_gap = n-1
std::vector<double> linspace (double min, double max, int n)
{
    // both min and max are stored in the vector
    if (n < 2)
    {
        std::cerr << "myMATH::Linspace: Cannot generate line space: "
                  << "n = " << n << ", "
                  << "less than 2!"
                  << std::endl;
        throw error_size;
    }

    double delta = (max - min) / (n - 1);

    std::vector<double> res(n);
    for (int i = 0; i < n; i ++)
    {
        res[i] = (min + delta * i);
    }

    return res;
};

// Trilinear interpolation
double triLinearInterp(AxisLimit const& grid_limit, std::vector<double> const& value, std::vector<double> const& pt_vec)
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
    
    // Key requirement:
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
        throw error_range;
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

// Create unit vector
Pt3D createUnitVector (Pt3D& pt1, Pt3D& pt2)
{
    Pt3D res = pt2 - pt1;
    res /= res.norm();
    return res;
};

// Create unit vector
Pt2D createUnitVector (Pt2D& pt1, Pt2D& pt2)
{
    Pt2D res = pt2 - pt1;
    res /= res.norm();
    return res;
};

// Calculate dot product
double dot (Pt3D const& pt1, Pt3D const& pt2)
{
    double res = 0;
    for (int i = 0; i < 3; i ++)
    {
        res += pt1[i] * pt2[i];
    }
    return res;
};

// Calculate dot product
double dot (Pt2D const& pt1, Pt2D const& pt2)
{
    double res = 0;
    for (int i = 0; i < 2; i ++)
    {
        res += pt1[i] * pt2[i];
    }
    return res;
};

// Calculate the distance between two points
double distance (Pt3D& pt1, Pt3D& pt2)
{
    Pt3D res = pt2 - pt1;
    return res.norm();
};

// Calculate the distance between two points
double distance (Pt2D& pt1, Pt2D& pt2)
{
    Pt2D res = pt2 - pt1;
    return res.norm();
};

// Calculate the distance between point and line
double distance (Pt3D& pt, Line3D& line)
{
    Pt3D diff = pt - line._pt;

    double dist_proj = dot(diff, line._unit_vector);
    double dist = std::pow(diff.norm(), 2) - dist_proj * dist_proj;
    
    if (dist >= 0)
    {
        dist = std::sqrt(dist);
    }
    else if (dist < 0 && dist > -SMALLNUMBER)
    {
        dist = 0;
    }
    else
    {
        std::cerr << "myMATH::Distance: negative distance" << std::endl;
        throw error_range;
    }

    return dist;
};

// Calculate the distance between point and line
double distance (Pt2D& pt, Line2D& line)
{
    Pt2D diff = pt - line._pt;

    double dist_proj = dot(diff, line._unit_vector);
    double dist = std::pow(diff.norm(), 2) - dist_proj * dist_proj;
    
    if (dist >= 0)
    {
        dist = std::sqrt(dist);
    }
    else if (dist < 0 && dist > -SMALLNUMBER)
    {
        dist = 0;
    }
    else
    {
        std::cerr << "myMATH::Distance: negative distance" << std::endl;
        throw error_range;
    }

    return dist;
};

}
