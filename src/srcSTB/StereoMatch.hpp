#ifndef STEREOMATCH_HPP
#define STEREOMATCH_HPP

#include "StereoMatch.h"


//---------------------//
// Auxiliary functions //
//---------------------//

// Find cross point of two lines 
//
//
template<class T>
Matrix<double> StereoMatch<T>::FindCrossPoint (Matrix<double> const& pt_1, Matrix<double> const& line_of_sight_1, Matrix<double> const& pt_2, Matrix<double> const& line_of_sight_2)
{
    // pt_center_1oncur + line_of_sight_1 * lambda_1 
    //  = pt_center_2oncur + line_of_sight_2 * lambda_2
    // pt_center_1oncur = (x1,x2,0); pt_center_2oncur = (x3,x4,0)
    // line_of_sight_1  = (c1,c2,0); line_of_sight_2  = (c3,c4,0)
    // lambda_1 = [(x3.c4-x4.c3) - (x1.c4-x2.c3)] / (c1.c4 - c2.c3)
    // lambda_2 = [(x2.c1-x1.c2) - (x4.c1-x3.c2)] / (c1.c4 - c2.c3)
    
    double num_0 = line_of_sight_1[0] * line_of_sight_2[1] - line_of_sight_1[1] * line_of_sight_2[0]; // c1.c4 - c2.c3

    // line_of_sight_1 is parallel to line_of_sight_2
    if (std::fabs(num_0) < MAGSMALLNUMBER)
    {
        std::cerr << "The two lines are parallel! "
                  << "No crossing point"
                  << std::endl;
        throw;
    }

    // line_of_sight_1 is not parallel to line_of_sight_2
    double num_1 = pt_2[0] * line_of_sight_2[1] - pt_2[1] * line_of_sight_2[0]; // x3.c4-x4.c3
    double num_2 = pt_1[0] * line_of_sight_2[1] - pt_1[1] * line_of_sight_2[0]; // x1.c4-x2.c3
    double lambda_1 = (num_1 - num_2) / num_0;
    // Matrix<double> pt_cross = pt_1 + line_of_sight_1 * lambda_1;
    Matrix<double> pt_cross (line_of_sight_1);
    pt_cross *= lambda_1;
    pt_cross += pt_1;

    return pt_cross;
}


// Given Y_pixel find X_pixel (col id)
// (Xu,Yu) is on the projected line 
//
template<class T>
int StereoMatch<T>::ComputeXpixel (Camera& cam, Matrix<double>& pt, Matrix<double>& line_of_sight, double y_pixel)
{
    // Tsai's model: xf,yf start from 1
    // But here, xf,yf start from 0
    // Note that y_pixel cannot be ref!
    y_pixel += 1;

    double y_d = cam.GetHpix() * (- y_pixel + cam.GetNpixh()/2.0 - cam.GetNoffh());
    double nx = line_of_sight[0];
    double ny = line_of_sight[1];
    double kr = cam.GetKr();
    double x_c = pt[0];
    double y_c = pt[1];
    int x_pixel;

    double x_d, x_d_1, x_d_part, x_d_2_num, x_d_2_den, x_d_3;

    if (kr > SMALLNUMBER)
    {
        x_d_1 = 2 * nx * y_d;

        x_d_part = std::sqrt( 
                std::pow(kr,3) 
                * ( - 4 * std::pow( 
                        kr * nx*nx * y_d*y_d 
                        - 3 * ny*ny * (1 + kr*y_d*y_d), 
                        3 ) 
                    + kr * std::pow(
                        27 * std::pow(ny,3) * x_c 
                        + 2 * kr * std::pow(nx,3) * std::pow(y_d,3) 
                        - 9 * nx * ny*ny 
                            * (3*y_c - 2*(y_d + kr * std::pow(y_d,3))), 
                        2 ) 
                  ) 
                );

        x_d_2_num = 2 * std::pow(2,1./3.) * (-kr * nx*nx * y_d*y_d + 3 * ny*ny * (1 + kr * y_d*y_d));
        x_d_2_den = std::pow( 
            - 27 * kr*kr * std::pow(ny,3) * x_c 
            + 27 * kr*kr * nx * ny*ny * y_c 
            - 18 * kr*kr * nx * ny*ny * y_d 
            - 2  * std::pow(kr,3) * std::pow(nx,3) * std::pow(y_d,3) 
            - 18 * std::pow(kr,3) * nx * ny*ny * std::pow(y_d,3) 
            + x_d_part, 
            1./3. );
        
        x_d_3 = std::pow(2, 2./3.) * x_d_2_den;

        x_d = 1 / (6*ny) * (x_d_1 + x_d_2_num/x_d_2_den - x_d_3 / kr);
    }
    else 
    {
        // Xd = Xu = Xc + nx/ny * (Yd-Yc)
        if (std::fabs(ny) > SMALLNUMBER)
        {
            x_d = x_c + nx / ny * (y_d - y_c);
        }
        else
        {
            std::cerr << "StereoMatch.ComputeXpixel: " 
                      << "ny of line of sight is zero, "
                      << "could not find unique x"
                      << std::endl;
            throw;
        }
    }

    // col index
    x_pixel = x_d/cam.GetWpix()/cam.GetKx() + cam.GetNpixw()/2 + cam.GetNoffw();

    // Tsai's model: xf,yf start from 1
    // But here, xf,yf start from 0
    x_pixel -= 1;

    return x_pixel;
}


// Given X_pixel find Y_pixel (row id)
// (Xu,Yu) is on the projected line 
//
template<class T>
int StereoMatch<T>::ComputeYpixel (Camera& cam, Matrix<double>& pt, Matrix<double>& line_of_sight, double x_pixel)
{
    // Tsai's model: xf,yf start from 1
    // But here, xf,yf start from 0
    x_pixel += 1;

    double x_d = x_pixel - cam.GetNpixw()/2.0 - cam.GetNoffw();
    x_d *= cam.GetWpix();
    x_d *= cam.GetKx();
    double nx = line_of_sight[0];
    double ny = line_of_sight[1];
    double kr = cam.GetKr();
    double x_c = pt[0];
    double y_c = pt[1];
    int y_pixel;

    double y_d, y_d_1, y_d_part, y_d_2_num, y_d_2_den, y_d_3;

    if (kr > SMALLNUMBER)
    {
        y_d_1 = 2 * ny * x_d;

        y_d_part = std::sqrt( 
                std::pow(kr,3) 
                * ( 4 * std::pow( 
                        - kr * ny*ny * x_d*x_d 
                        + 3 * nx*nx * (1 + kr*x_d*x_d), 
                        3 ) 
                    + kr * std::pow(
                        27 * std::pow(nx,3) * y_c 
                        + 2 * kr * std::pow(ny,3) * std::pow(x_d,3) 
                        - 9 * nx*nx * ny 
                            * (3*x_c - 2*(x_d + kr * std::pow(x_d,3))), 
                        2 ) 
                  ) 
                );

        y_d_2_num = 2 * std::pow(2,1./3.) * (-kr * ny*ny * x_d*x_d + 3 * nx*nx * (1 + kr * x_d*x_d));
        y_d_2_den = std::pow( 
            - 27 * kr*kr * nx*nx * ny * x_c 
            + 27 * kr*kr * std::pow(nx,3) * y_c 
            + 18 * kr*kr * nx*nx * ny * x_d 
            + 2  * std::pow(kr,3) * std::pow(ny,3) * std::pow(x_d,3) 
            + 18 * std::pow(kr,3) * nx*nx * ny * std::pow(x_d,3) 
            + y_d_part, 
            1./3. );
        
        y_d_3 = std::pow(2,2./3.) * y_d_2_den;

        y_d = 1 / (6*nx) * (y_d_1 - y_d_2_num/y_d_2_den + y_d_3/kr);
    }
    else 
    {
        // Yd = Yu = Yc + ny/nx * (Xd-Xc)
        if (std::fabs(nx) > SMALLNUMBER)
        {
            y_d = y_c + ny / nx * (x_d - x_c);
        }
        else
        {
            std::cerr << "StereoMatch.ComputeXpixel: " 
                      << "nx of line of sight is zero, "
                      << "could not find unique y"
                      << std::endl;
            throw;
        }
    }

    y_pixel = - y_d/cam.GetHpix() + cam.GetNpixh()/2 - cam.GetNoffh();

    // Tsai's model: xf,yf start from 1
    // But here, xf,yf start from 0
    y_pixel -= 1;

    return y_pixel;
}


template<class T>
void StereoMatch<T>::Triangulation(int n_cam_tri, std::vector<Matrix<double>>& pt_2d_list, Matrix<double>& pt_world, double& error)
{
    Matrix<double> mtx(3, 3, 0);
    Matrix<double> pt_3d(3, 1, 0);

    std::vector<Matrix<double>> line_of_sight_list;

    for (int i = 0; i < n_cam_tri; i++)
    {
        // vector from camera center to image pt in world coordinate
        Matrix<double> line_of_sight(
            _cam_list[i].ImgMMToWorld(pt_2d_list[i]),
            _cam_list[i].GetCenterInWorld()
        );
        line_of_sight_list.push_back(line_of_sight);

        double nx = line_of_sight[0];
        double ny = line_of_sight[1];
        double nz = line_of_sight[2];

        Matrix<double> tmp(
            1 - nx * nx, -nx * ny, -nx * nz,
            -nx * ny, 1 - ny * ny, -ny * nz,
            -nz * nx, -nz * ny, 1 - nz * nz
        );

        pt_3d += tmp * _cam_list[i].GetCenterInWorld();
        mtx += tmp;
    }

    pt_world = mtx.Inverse() * pt_3d;

    // calculate errors 
    error = 0;
    Matrix<double> diff_vec(3,1);
    for (int i = 0; i < n_cam_tri; i++)
    {            
        error += pt_world.Dist(_cam_list[i].GetCenterInWorld(), line_of_sight_list[i]);
    }
    error /= n_cam_tri;
};


//
// Major functions
//

template<class T>
StereoMatch<T>::StereoMatch(
    std::vector<Camera> const& cam_list,
    double tor_2d,
    double tor_3d
) : _cam_list(cam_list), _n_cam(cam_list.size()), _tor_2d(tor_2d), _tor_3d(tor_3d)
{}


template<class T>
void StereoMatch<T>::SetObjectListMMFromMM (std::vector<std::vector<T>>& object_list_mm)
{
    _object_list_mm = object_list_mm;
}


template<class T>
void StereoMatch<T>::SetObjectListMMFromPixel (std::vector<std::vector<T>>& object_list_pixel)
{
    clock_t t_start, t_end;

    t_start = clock();

    _object_list_mm = object_list_pixel;
    Matrix<double> pt_mm(3,1);

    for (int cam_id = 0; cam_id < _object_list_mm.size(); cam_id ++)
    {
        for (int object_id = 0; object_id < _object_list_mm[cam_id].size(); object_id ++)
        {
            pt_mm = _object_list_mm[cam_id][object_id].GetCenterPos();
            pt_mm = _cam_list[cam_id].ImgPixelToMM(pt_mm);
            _object_list_mm[cam_id][object_id].SetCenterPos(pt_mm);
        }
    }

    t_end = clock();
    std::cout << "SetObjectListMMFromPixel: " 
              << (double) (t_end - t_start)/CLOCKS_PER_SEC
              << "s"
              << std::endl;

}


template<class T>
void StereoMatch<T>::ClearAll ()
{
    ClearCamList();     
    ClearObjectListMM();     
    ClearObjectIDMap();     
    ClearObjectIDMatchList();
}


template<class T>
void StereoMatch<T>::MakeObjectIDMap(std::vector<std::vector<T>>& object_list_pixel)
{
    clock_t t_start, t_end;

    t_start = clock();

    for (int i = 0; i < _n_cam; i++)
    {
        // note that:
        // _n_pix_w: number of pixel along X direction (col)
        // _n_pix_h: number of pixel along Y direction (row)
        int n_row = _cam_list[i].GetNpixh();
        int n_col = _cam_list[i].GetNpixw();

        #pragma region InitializationMap
        // initialize for object_id_map_one (each camera)
        // n_row * n_col * object count
        std::vector<std::vector<std::vector<int>>> object_id_map_one;

        if (i == 0)
        {
            // no need to create a map for 1st cam 
            _object_id_map.push_back(object_id_map_one);
            continue;
        }
        // initialize for row_object_id 
        // n_col * object count
        object_id_map_one.reserve(n_row+2);

        std::vector<std::vector<int>> row_object_id;
        row_object_id.reserve(n_col+2);

        std::vector<int> grid_object_id;
        grid_object_id.reserve(20);
        grid_object_id.push_back(-1);
        for (int col_id = 0; col_id < n_col; col_id ++)
        {
            row_object_id.push_back(grid_object_id);
        }

        for (int row_id = 0; row_id < n_row; row_id++)
        {
            object_id_map_one.push_back(row_object_id);
        }
        #pragma endregion initialize _object_id_map
        
        #pragma region FillMap
        for (int j = 0; j < object_list_pixel[i].size(); j++)
        {
            int row = object_list_pixel[i][j].GetCenterPos()[1]; // y coordinate
            int col = object_list_pixel[i][j].GetCenterPos()[0]; // x coordinate 

            // check whether (x,y)=(col,row) is within the domain 
            if (row > n_row - 1 || col > n_col - 1 || row < 0 || col < 0)
            {
                continue;
            }

            // check whether there is already a object stored in this pixel
            if (object_id_map_one[row][col][0] != -1) 
            {
                object_id_map_one[row][col].push_back(j);
            }
            else
            {
                object_id_map_one[row][col][0] = j;
            }
        }
        _object_id_map.push_back(object_id_map_one);
        #pragma endregion
    }
    
    t_end = clock();
    std::cout << "MakeObjectIDMap: " 
              << (double) (t_end - t_start)/CLOCKS_PER_SEC
              << "s"
              << std::endl;
}

// TBC
template<class T>
PixelRange StereoMatch<T>::FindSearchRegion (int cam_id, std::vector<Matrix<double>>& pt_list, std::vector<Matrix<double>>& unit_list)
{
    Camera cam = _cam_list[cam_id];
    PixelRange search_region;
    int search_region_state = 0;

    if (pt_list.size() == 1)
    {
        search_region._row_max = cam.GetNpixh();
        search_region._col_max = cam.GetNpixw();
        return search_region;
    }

    // Find all crossing points
    Matrix<double> pt_cur(3,1);
    Matrix<double> unit_cur(3,1);
    Matrix<double> pt_test(3,1);
    Matrix<double> unit_test(3,1);
    Matrix<double> per_unit_cur(3,1);
    Matrix<double> per_unit_test(3,1);
    Matrix<double> line_cur_plus(3,1);
    Matrix<double> line_cur_minus(3,1);
    Matrix<double> line_test_plus(3,1);
    Matrix<double> line_test_minus(3,1);
    Matrix<double> pt_cross(3,1);

    for (int i = 0; i < int(pt_list.size())-1; i ++)
    {
        pt_cur = pt_list[i];
        unit_cur = unit_list[i];
        per_unit_cur(0,0) = - unit_cur[1];
        per_unit_cur(1,0) = unit_cur[0];
        per_unit_cur(2,0) = 0;

        for (int j = i+1; j < pt_list.size(); j ++)
        {
            pt_test = pt_list[j];
            unit_test = unit_list[j];
            per_unit_test(0,0) = - unit_test[1];
            per_unit_test(1,0) = unit_test[0];
            per_unit_test(2,0) = 0;

            line_cur_plus = pt_cur + per_unit_cur * _tor_2d;
            line_cur_minus = pt_cur - per_unit_cur * _tor_2d;
            line_test_plus = pt_test + per_unit_test * _tor_2d;
            line_test_minus = pt_test - per_unit_test * _tor_2d;

            pt_cross = FindCrossPoint(
                line_cur_plus, 
                unit_cur,
                line_test_plus,
                unit_test
            );
            pt_cross = cam.ImgMMToPixel(pt_cross);
            if (search_region_state == 0)
            {
                search_region._row_min = int(pt_cross[1]);
                search_region._row_max = search_region._row_min + 1;
                search_region._col_min = int(pt_cross[0]);
                search_region._col_max = search_region._col_min + 1;

                search_region_state = 1;
            }
            else 
            {
                search_region.SetRange(pt_cross[1], pt_cross[0]);
            }
            

            pt_cross = FindCrossPoint(
                line_cur_minus, 
                unit_cur,
                line_test_plus,
                unit_test
            );
            pt_cross = cam.ImgMMToPixel(pt_cross);
            search_region.SetRange(pt_cross[1], pt_cross[0]);

            pt_cross = FindCrossPoint(
                line_cur_plus, 
                unit_cur,
                line_test_minus,
                unit_test
            );
            pt_cross = cam.ImgMMToPixel(pt_cross);
            search_region.SetRange(pt_cross[1], pt_cross[0]);

            pt_cross = FindCrossPoint(
                line_cur_minus, 
                unit_cur,
                line_test_minus,
                unit_test
            );
            pt_cross = cam.ImgMMToPixel(pt_cross);
            search_region.SetRange(pt_cross[1], pt_cross[0]);

        }
    }

    return search_region;
}


template<class T>
PixelRange StereoMatch<T>::FindSearchRegionOrig (int cam_id, std::vector<Matrix<double>>& pt_list, std::vector<Matrix<double>>& unit_list)
{
    Camera cam = _cam_list[cam_id];
    PixelRange search_region;

    if (pt_list.size() == 1)
    {
        // cam_id = 1, 2nd cam
        search_region._row_max = cam.GetNpixh();
        search_region._col_max = cam.GetNpixw();
        return search_region;
    }
    else if (pt_list.size() == 2) 
    {
        // cam_id = 2, 3nd cam
        Matrix<double> pt_cross (
            FindCrossPoint (
                pt_list[0], 
                unit_list[0], 
                pt_list[1], 
                unit_list[1]
            )
        );

        double cos = std::max(std::sqrt(std::fabs(unit_list[0].Dot(unit_list[1]))), MAGSMALLNUMBER);
        double sin = std::max(std::sqrt(1-cos*cos), MAGSMALLNUMBER);

        double width = 1.0;
        if (cos < sin)
        {
            width = _tor_2d/cos;
        }
        else 
        {
            width = _tor_2d/sin;
        }
        width = std::min(width, 1.2*_fov_length_mm);

        Matrix<double> pt_upright(pt_cross);
        pt_upright(0,0) = pt_upright(0,0) + width;
        pt_upright(1,0) = pt_upright(1,0) + width;
        pt_upright = cam.ImgMMToPixel(pt_upright);

        Matrix<double> pt_lowleft(pt_cross);
        pt_lowleft(0,0) = pt_lowleft(0,0) - width;
        pt_lowleft(1,0) = pt_lowleft(1,0) - width;
        pt_lowleft = cam.ImgMMToPixel(pt_lowleft);

        search_region._row_min = std::floor(pt_upright(1,0));
        search_region._row_max = std::ceil(pt_lowleft(1,0))+1;
        search_region._col_min = std::floor(pt_lowleft(0,0));
        search_region._col_max = std::ceil(pt_upright(0,0))+1;

    }

    return search_region;
}


template<class T>
PixelRange StereoMatch<T>::FindSearchRegion (int cam_id, int line_ref_id, std::vector<Matrix<double>>& pt_list, std::vector<Matrix<double>>& unit_list)
{
    if (line_ref_id > int(pt_list.size())-1)
    {
        std::cout << "FindSearchRegion: "
                  << "line reference id is: "
                  << line_ref_id << ", "
                  << "out of range: "
                  << pt_list.size()
                  << std::endl;
        throw;
    }

    Camera cam = _cam_list[cam_id];
    PixelRange search_region;
    int search_region_state = 0;

    if (pt_list.size() == 1)
    {
        search_region._row_max = cam.GetNpixh();
        search_region._col_max = cam.GetNpixw();
        return search_region;
    }

    // Find all crossing points
    Matrix<double> pt_cur (pt_list[line_ref_id]);
    Matrix<double> unit_cur (unit_list[line_ref_id]);
    Matrix<double> per_unit_cur (
        3,1,
        - unit_cur[1],
        unit_cur[0],
        0
    );
    Matrix<double> pt_test(3,1);
    Matrix<double> unit_test(3,1);
    Matrix<double> per_unit_test(3,1);
    Matrix<double> pt_cross(3,1);

    for (int j = 0; j < pt_list.size(); j ++)
    {
        if (j == line_ref_id)
        {
            continue;
        }

        pt_test = pt_list[j];
        unit_test = unit_list[j];
        per_unit_test(0,0) = - unit_test[1];
        per_unit_test(1,0) = unit_test[0];
        per_unit_test(2,0) = 0;

        pt_cross = FindCrossPoint(
            pt_cur + per_unit_cur * _tor_2d, 
            unit_cur,
            pt_test + per_unit_test * _tor_2d,
            unit_test
        );
        pt_cross = cam.ImgMMToPixel(pt_cross);
        if (search_region_state == 0)
        {
            search_region._row_min = int(pt_cross[1]);
            search_region._row_max = search_region._row_min + 1;
            search_region._col_min = int(pt_cross[0]);
            search_region._col_max = search_region._col_min + 1;

            search_region_state = 1;
        }
        else 
        {
            search_region.SetRange(pt_cross[1], pt_cross[0]);
        }

        pt_cross = FindCrossPoint(
            pt_cur - per_unit_cur * _tor_2d, 
            unit_cur,
            pt_test + per_unit_test * _tor_2d,
            unit_test
        );
        pt_cross = cam.ImgMMToPixel(pt_cross);
        search_region.SetRange(pt_cross[1], pt_cross[0]);

        pt_cross = FindCrossPoint(
            pt_cur + per_unit_cur * _tor_2d, 
            unit_cur,
            pt_test - per_unit_test * _tor_2d,
            unit_test
        );
        pt_cross = cam.ImgMMToPixel(pt_cross);
        search_region.SetRange(pt_cross[1], pt_cross[0]);

        pt_cross = FindCrossPoint(
            pt_cur - per_unit_cur * _tor_2d, 
            unit_cur,
            pt_test - per_unit_test * _tor_2d,
            unit_test
        );
        pt_cross = cam.ImgMMToPixel(pt_cross);
        search_region.SetRange(pt_cross[1], pt_cross[0]);

    }
    

    return search_region;
}


template<class T>
void StereoMatch<T>::ObjectMarkerAroundLine(int cam_id, Matrix<double>& pt, Matrix<double>& line_of_sight, Matrix<int>& object_marker, PixelRange& search_region)
{
    Camera cam = _cam_list[cam_id];
    Matrix<double> per_line_of_sight (
        3,1,
        -line_of_sight[1], line_of_sight[0], 0
    );

    int n_pix_h = cam.GetNpixh(); // number of rows
    int n_pix_w = cam.GetNpixw(); // number of cols

    // Search the candidate within the torlerance near a line
    // Determine the border of the torlerance by iterate each x pixel or y pixel depending on the slope.
    // x_pixel (col id), y_pixel (row id)
    int x_pixel_1, x_pixel_2, y_pixel_1, y_pixel_2; 
    int min, max;

    // Matrix<double> pt_plus  = pt + per_line_of_sight*_tor_2d;
    // Matrix<double> pt_minus = pt - per_line_of_sight*_tor_2d;
    Matrix<double> pt_plus (
        3,1,
        pt[0] + per_line_of_sight[0]*_tor_2d,
        pt[1] + per_line_of_sight[1]*_tor_2d,
        0
    );
    Matrix<double> pt_minus (
        3,1,
        pt[0] - per_line_of_sight[0]*_tor_2d,
        pt[1] - per_line_of_sight[1]*_tor_2d,
        0
    );

    // determine the search orientation on the pixel image 
    //  if the |slope| >= 1, iterate every y_pixel (row)
    //  else, iterate every x_pixel (col) 
    if (std::fabs(line_of_sight[0]) == 0 || std::fabs(line_of_sight[1]/line_of_sight[0]) >= 1)
    {
        // iterate every y_pixel (row id)
        // Xu = Xd (1 + k1 Xd^2 + k1 Yd^2)
        // Yu = Yd (1 + k1 Xd^2 + k1 Yd^2)
        // (Xu  = (Xc  + lambda (nx   
        //  Yu)    Yc)           ny)
        // using mathematica to compute the results     
        int y_start = std::max(0, std::min(search_region._row_min, n_pix_h-1));
        int y_end   = std::min(n_pix_h, std::max(search_region._row_max, 1));
        for (int y_pixel = y_start; y_pixel < y_end; y_pixel ++)
        {
            // Get x_pixel (col id) from two lines 
            x_pixel_1 = ComputeXpixel (
                cam, 
                pt_plus, 
                line_of_sight, 
                y_pixel 
            );
            x_pixel_2 = ComputeXpixel (
                cam, 
                pt_minus, 
                line_of_sight, 
                y_pixel 
            );
            if (x_pixel_1 > x_pixel_2)
            {
                max = x_pixel_1 + 1; 
                min = x_pixel_2;
            }
            else 
            {
                max = x_pixel_2 + 1; 
                min = x_pixel_1;
            }
            
            min = std::max(0, min);
            min = std::min(n_pix_w-1, min);
            max = std::max(1, max);
            max = std::min(n_pix_w, max);

            for (int col = min; col < max; col ++)
            {
                if (_object_id_map[cam_id][y_pixel][col][0] != -1) 
                {
                    object_marker(y_pixel, col) += 1;
                }
            }
        }
    }
    else 
    {
        // iterate every x_pixel (col id)
        // Xu = Xd (1 + k1 Xd^2 + k1 Yd^2)
        // Yu = Yd (1 + k1 Xd^2 + k1 Yd^2)
        // (Xu  = (Xc  + lambda (nx   
        //  Yu)    Yc)           ny)
        // using mathematica to compute the results
        int x_start = std::max(0, std::min(search_region._col_min, n_pix_w-1));
        int x_end   = std::min(n_pix_w, std::max(search_region._col_max, 1));
        for (int x_pixel = x_start; x_pixel < x_end; x_pixel ++)
        {
            // Get x_pixel (col id) from two lines 
            y_pixel_1 = ComputeYpixel (
                cam, 
                pt_plus, 
                line_of_sight, 
                x_pixel 
            );
            y_pixel_2 = ComputeYpixel (
                cam, 
                pt_minus, 
                line_of_sight, 
                x_pixel 
            );
            if (y_pixel_1 > y_pixel_2)
            {
                max = y_pixel_1 + 1; 
                min = y_pixel_2;
            }
            else 
            {
                max = y_pixel_2 + 1; 
                min = y_pixel_1;
            }

            min = std::max(0, min);
            min = std::min(n_pix_h-1, min);
            max = std::max(0, max);
            max = std::min(n_pix_h, max);

            for (int row = min; row < max; row ++)
            {
                if (_object_id_map[cam_id][row][x_pixel][0] != -1) 
                {
                    object_marker(row, x_pixel) += 1;
                }
            }
        }
    }
}


template<class T>
void StereoMatch<T>::FindTracerMatch (
    int cam_cur_id, 
    std::deque<std::vector<int>>& tracer_id_match_list, 
    std::vector<int> tracer_id_match, 
    std::deque<double>& error_list
)
{
    int n_row = _cam_list[cam_cur_id].GetNpixh();
    int n_col = _cam_list[cam_cur_id].GetNpixw();

    if (cam_cur_id < 1)
    {
        std::cerr << "StereoMatch error: "
                  << "camera id is too small to "
                  << "start finding tracer match!"
                  << std::endl;
        throw;
    }
    
    // tracer marker around the line of sight 
    // intialize a 2D map with value 0
    Matrix<int> object_marker(
        n_row, n_col,
        0 
    );

    std::vector<Matrix<double>> line_of_sight_list;
    std::vector<Matrix<double>> pt_list;

    // Mark the intersection from all previous cameras on the current camera
    #pragma region FindSearchRegion
    for (int i = 0; i < cam_cur_id; i ++)
    {       
        // project from cam_i onto cam_cur
        // find out tracer position in camera coordinate [mm]
        Matrix<double> pt_mm = _object_list_mm[i][tracer_id_match[i]].GetCenterPos();
        // find out tracer position in world coordinate
        Matrix<double> pt_world = _cam_list[i].ImgMMToWorld(pt_mm);

        // position of cam_i center on cam_cur
        Matrix<double> pt_center_i_on_cur = _cam_list[cam_cur_id].WorldToImgMM(_cam_list[i].GetCenterInWorld());
        pt_list.push_back(pt_center_i_on_cur);
        // position of tracer center on cam_i projected on cam_cur
        Matrix<double> pt_tracer_i_on_cur = _cam_list[cam_cur_id].WorldToImgMM(pt_world);

        // unit vector in projected line of sight direction 
        Matrix<double> line_of_sight(pt_tracer_i_on_cur - pt_center_i_on_cur);
        line_of_sight /= line_of_sight.Mag();
        line_of_sight_list.push_back(line_of_sight);
    }
    #pragma endregion Find Search Region

    #pragma region MarkProjection
    for (int i = 0; i < cam_cur_id; i ++)
    {
        PixelRange search_region = FindSearchRegion(cam_cur_id, i, pt_list, line_of_sight_list);

        // find tracer candidates around the line of sight 
        // increase object_marker by 1 if found 
        ObjectMarkerAroundLine(
            cam_cur_id,
            pt_list[i], 
            line_of_sight_list[i],
            object_marker,
            search_region
        );   
    }
    #pragma endregion project Line of sight and mark tracers near the line

    PixelRange search_region_all = FindSearchRegion(cam_cur_id, pt_list, line_of_sight_list);
    int row_start = std::min(std::max(search_region_all._row_min, 0), n_row-1);
    int row_end   = std::min(std::max(search_region_all._row_max, 1), n_row);
    int col_start = std::min(std::max(search_region_all._col_min, 0), n_col-1);
    int col_end   = std::min(std::max(search_region_all._col_max, 1), n_col);
    // int row_start = 0;
    // int row_end   = n_row;
    // int col_start = 0;
    // int col_end   = n_col;

    for (int i = row_start; i < row_end; i ++)
    {
        for (int j = col_start; j < col_end; j ++)
        {
            if (object_marker(i,j) == cam_cur_id)
            {
                // judge whether the distances between the candidate
                // and all the lines are all within the range 
                for (
                    int k = 0; 
                    k < _object_id_map[cam_cur_id][i][j].size(); 
                    k ++
                )
                {
                    int tracer_id;
                    tracer_id = _object_id_map[cam_cur_id][i][j][k];

                    if (tracer_id == -1)
                    {
                        std::cerr << "tracer_id is wrong" << std::endl;
                        std::cerr << "(cam_cur_id,i,j,k) = " << "(" 
                                  << cam_cur_id << ", " 
                                  << i << ", "
                                  << j << ", "
                                  << k << ")"
                                  << std::endl;
                        throw;
                    }

                    tracer_id_match.push_back(tracer_id);

                    #pragma region CheckTriangulation
                    // test 3d distance
                    // by triangulation                 
                    std::vector<Matrix<double>> pt_2d_list;
                    for (int m = 0; m < cam_cur_id+1; m ++)
                    {
                        Matrix<double> pt_2d(_object_list_mm[m][tracer_id_match[m]].GetCenterPos());
                        pt_2d_list.push_back(pt_2d);
                    }

                    Matrix<double> pt_world;
                    double error;

                    Triangulation(
                        cam_cur_id+1,
                        pt_2d_list,
                        pt_world,
                        error
                    );
                    #pragma endregion

                    // check 3D here, triangulation
                    if (error > _tor_3d) 
                    {
                        tracer_id_match.pop_back(); // if not passing the error, then delete the candidate
                    } 
                    else 
                    {
                        // if there is still other cameras to check
                        if (cam_cur_id >= 1 && cam_cur_id < _check_cam_id-1) // the first condition is unnecessary?
                        {
                            int cam_next_id = cam_cur_id + 1;

                            // move onto the next camera and search candidates
                            FindTracerMatch (
                                cam_next_id,
                                tracer_id_match_list,
                                tracer_id_match, 
                                error_list                  
                            );
                            tracer_id_match.pop_back(); // clear the match for the next candidate.
                        }
                        // if the current camera is the last one, then finalize the match.
                        else
                        {
                            if (_check_cam_id < _n_cam)
                            {
                                int cam_next_id = cam_cur_id + 1;

                                CheckTracerMatch (
                                    cam_next_id,
                                    tracer_id_match_list,
                                    tracer_id_match, 
                                    pt_world,
                                    error_list
                                );

                                tracer_id_match.pop_back();
                            }
                            else 
                            {
                                // // Print info 
                                // std::cout << "Newly found match list: ";
                                // for (int id = 0; id < _n_cam; id ++)
                                // {
                                //     std::cout << tracer_id_match[id] 
                                //               << " ";
                                // }
                                // std::cout << std::endl;

                                // tracer_id_match_list.push_back(std::move(tracer_id_match));
                                tracer_id_match_list.push_back(tracer_id_match);

                                error_list.push_back(error);

                                tracer_id_match.pop_back();
                            }
                        }
                    }
                }
            }
        }
    }
}


// pending
template<class T>
void StereoMatch<T>::FindTracerMatchOrig (
    int cam_cur_id,
    std::deque<std::vector<int>>& tracer_id_match_list,  
    std::vector<int> tracer_id_match,
    std::deque<double>& error_list
)
{
    int n_row = _cam_list[cam_cur_id].GetNpixh();
    int n_col = _cam_list[cam_cur_id].GetNpixw();

    //                          //
    // cam_cur_id = 1 (2nd cam) //
    //                          //
    // tracer marker around the line of sight 
    // intialize a 2D map with value 0
    // optimization for cam_cur_id == 1
    if (cam_cur_id == 1)  
    {
        #pragma region MarkProjection

        // find tracer candidates around the line of sight         
        // ObjectMarkerAroundLine   
        Camera cam = _cam_list[cam_cur_id];
        Matrix<double> pt (
            _cam_list[cam_cur_id].WorldToImgMM(_cam_list[0].GetCenterInWorld())
        ); 
        Matrix<double> line_of_sight(
            _cam_list[cam_cur_id].WorldToImgMM(_cam_list[0].GetCenterInWorld()), 
            _cam_list[cam_cur_id].WorldToImgMM(
                _cam_list[0].ImgMMToWorld(
                    _object_list_mm[0][tracer_id_match[0]].GetCenterPos()
                )
            )
        );
        Matrix<double> per_line_of_sight (
            3,1,
            -line_of_sight[1], line_of_sight[0], 0
        );

        // Search the candidate within the torlerance near a line
        // Determine the border of the torlerance by iterate each x pixel or y pixel depending on the slope.
        // x_pixel (col id), y_pixel (row id)
        int x_pixel_1, x_pixel_2, y_pixel_1, y_pixel_2; 
        int min, max;

        // Matrix<double> pt_plus  = pt + per_line_of_sight*_tor_2d;
        // Matrix<double> pt_minus = pt - per_line_of_sight*_tor_2d;
        Matrix<double> pt_plus (
            3,1,
            pt[0] + per_line_of_sight[0]*_tor_2d,
            pt[1] + per_line_of_sight[1]*_tor_2d,
            0
        );
        Matrix<double> pt_minus (
            3,1,
            pt[0] - per_line_of_sight[0]*_tor_2d,
            pt[1] - per_line_of_sight[1]*_tor_2d,
            0
        );


        // determine the search orientation on the pixel image 
        //  if the |slope| >= 1, iterate every y_pixel (row)
        //  else, iterate every x_pixel (col) 
        if (std::fabs(line_of_sight[0]) == 0 || std::fabs(line_of_sight[1]/line_of_sight[0]) >= 1)
        {
            // iterate every y_pixel (row id)
            // Xu = Xd (1 + k1 Xd^2 + k1 Yd^2)
            // Yu = Yd (1 + k1 Xd^2 + k1 Yd^2)
            // (Xu  = (Xc  + lambda (nx   
            //  Yu)    Yc)           ny)
            // using mathematica to compute the results     
            int y_start = 0;
            int y_end   = n_row;
            for (int y_pixel = y_start; y_pixel < y_end; y_pixel ++)
            {
                // Get x_pixel (col id) from two lines 
                x_pixel_1 = ComputeXpixel (
                    cam, 
                    pt_plus, 
                    line_of_sight, 
                    y_pixel 
                );
                x_pixel_2 = ComputeXpixel (
                    cam, 
                    pt_minus, 
                    line_of_sight, 
                    y_pixel 
                );
                if (x_pixel_1 > x_pixel_2)
                {
                    max = x_pixel_1 + 1; 
                    min = x_pixel_2;
                }
                else 
                {
                    max = x_pixel_2 + 1; 
                    min = x_pixel_1;
                }
                
                min = std::max(0, min);
                min = std::min(n_col-1, min);
                max = std::max(1, max);
                max = std::min(n_col, max);

                for (int col = min; col < max; col ++)
                {
                    if (_object_id_map[cam_cur_id][y_pixel][col][0] != -1) 
                    {
                        int size = _object_id_map[cam_cur_id][y_pixel][col].size();
                        for (int k = 0; k < size; k ++)
                        {
                            int tracer_id;
                            tracer_id = _object_id_map[cam_cur_id][y_pixel][col][k];

                            if (tracer_id == -1)
                            {
                                std::cerr << "tracer_id is wrong" << std::endl;
                                std::cerr << "(cam_cur_id,i,j,k) = " << "(" 
                                        << cam_cur_id << ", " 
                                        << y_pixel << ", "
                                        << col << ", "
                                        << k << ")"
                                        << std::endl;
                                throw;
                            }
                            
                            bool in_range = true;
                            double dist = _object_list_mm[cam_cur_id][tracer_id].GetCenterPos().Dist(pt, line_of_sight);

                            if (dist > _tor_2d)
                            {
                                in_range = false;
                            }

                            if (in_range)
                            {
                                in_range = CheckReProject(
                                    cam_cur_id, tracer_id_match, 
                                    _object_list_mm[cam_cur_id][tracer_id].GetCenterPos()
                                );
                            }

                            if (in_range)
                            {
                                tracer_id_match.push_back(tracer_id);

                                // if there is still other cameras to check
                                if (cam_cur_id >= 1 && cam_cur_id < _check_cam_id - 1) 
                                {
                                    int cam_next_id = cam_cur_id + 1;

                                    // move onto the next camera and search candidates
                                    FindTracerMatchOrig (
                                        cam_next_id,
                                        tracer_id_match_list,
                                        tracer_id_match, 
                                        error_list                  
                                    );

                                    tracer_id_match.pop_back(); // clear the match for the next candidate.
                                }
                                // if the current camera is the last one, then finalize the match.
                                else
                                {
                                    // std::cout << "cam_cur_id=" << cam_cur_id << ", "
                                    //           << "wrong!\n";
                                    // throw;

                                    //test 3d distance by triangulation                 
                                    std::vector<Matrix<double>> pt_2d_list;
                                    for (int m = 0; m < cam_cur_id+1; m ++)
                                    {
                                        if (m == cam_cur_id)
                                        {
                                            pt_2d_list.push_back(
                                                _object_list_mm[m][tracer_id].GetCenterPos()
                                            ); 
                                        }
                                        else 
                                        {
                                            pt_2d_list.push_back(
                                                _object_list_mm[m][tracer_id_match[m]].GetCenterPos()
                                            );
                                        } 
                                    }

                                    Matrix<double> pt_world(3,1);
                                    double error_3d=0.0;

                                    Triangulation(
                                        cam_cur_id+1,
                                        pt_2d_list,
                                        pt_world,
                                        error_3d
                                    );

                                    if (error_3d > _tor_3d)
                                    {
                                        tracer_id_match.pop_back();
                                    }
                                    else 
                                    {
                                        if (_check_cam_id < _n_cam)
                                        {
                                            int cam_next_id = cam_cur_id + 1;

                                            CheckTracerMatch (
                                                cam_next_id,
                                                tracer_id_match_list,
                                                tracer_id_match, 
                                                pt_world,
                                                error_list
                                            );

                                            tracer_id_match.pop_back();
                                        }
                                        else 
                                        {
                                            tracer_id_match_list.push_back(tracer_id_match);

                                            error_list.push_back(error_3d);

                                            tracer_id_match.pop_back();
                                        }
                                    }
                                    
                                }
                            }
                            
                        }
                    }
                }
            }
        }
        else 
        {
            // iterate every x_pixel (col id)
            // Xu = Xd (1 + k1 Xd^2 + k1 Yd^2)
            // Yu = Yd (1 + k1 Xd^2 + k1 Yd^2)
            // (Xu  = (Xc  + lambda (nx   
            //  Yu)    Yc)           ny)
            // using mathematica to compute the results
            int x_start = 0;
            int x_end   = n_col;
            for (int x_pixel = x_start; x_pixel < x_end; x_pixel ++)
            {
                // Get x_pixel (col id) from two lines 
                y_pixel_1 = ComputeYpixel (
                    cam, 
                    pt_plus, 
                    line_of_sight, 
                    x_pixel 
                );
                y_pixel_2 = ComputeYpixel (
                    cam, 
                    pt_minus, 
                    line_of_sight, 
                    x_pixel 
                );
                if (y_pixel_1 > y_pixel_2)
                {
                    max = y_pixel_1 + 1; 
                    min = y_pixel_2;
                }
                else 
                {
                    max = y_pixel_2 + 1; 
                    min = y_pixel_1;
                }

                min = std::max(0, min);
                min = std::min(n_row-1, min);
                max = std::max(0, max);
                max = std::min(n_row, max);

                for (int row = min; row < max; row ++)
                {
                    if (_object_id_map[cam_cur_id][row][x_pixel][0] != -1) 
                    {
                        int size = _object_id_map[cam_cur_id][row][x_pixel].size();
                        for (int k = 0; k < size; k ++)
                        {
                            int tracer_id;
                            tracer_id = _object_id_map[cam_cur_id][row][x_pixel][k];

                            if (tracer_id == -1)
                            {
                                std::cerr << "tracer_id is wrong" << std::endl;
                                std::cerr << "(cam_cur_id,i,j,k) = " << "(" 
                                        << cam_cur_id << ", " 
                                        << row << ", "
                                        << x_pixel << ", "
                                        << k << ")"
                                        << std::endl;
                                throw;
                            }

                            bool in_range = true;
                            double dist = _object_list_mm[cam_cur_id][tracer_id].GetCenterPos().Dist(pt, line_of_sight);

                            if (dist > _tor_2d)
                            {
                                in_range = false;
                            }

                            if (in_range)
                            {
                                in_range = CheckReProject(
                                    cam_cur_id, tracer_id_match, 
                                    _object_list_mm[cam_cur_id][tracer_id].GetCenterPos()
                                );
                            }

                            if (in_range)
                            {
                                tracer_id_match.push_back(tracer_id);

                                // if there is still other cameras to check
                                if (cam_cur_id >= 1 && cam_cur_id < _check_cam_id - 1) 
                                {
                                    int cam_next_id = cam_cur_id + 1;

                                    // move onto the next camera and search candidates
                                    FindTracerMatchOrig (
                                        cam_next_id,
                                        tracer_id_match_list,
                                        tracer_id_match, 
                                        error_list                  
                                    );

                                    tracer_id_match.pop_back(); // clear the match for the next candidate.
                                }
                                // if the current camera is the last one, then finalize the match.
                                else
                                {
                                    // std::cout << "cam_cur_id=" << cam_cur_id << ", "
                                    //           << "wrong!\n";
                                    // throw;

                                    //test 3d distance by triangulation                 
                                    std::vector<Matrix<double>> pt_2d_list;
                                    for (int m = 0; m < cam_cur_id+1; m ++)
                                    {
                                        if (m == cam_cur_id)
                                        {
                                            pt_2d_list.push_back(
                                                _object_list_mm[m][tracer_id].GetCenterPos()
                                            ); 
                                        }
                                        else 
                                        {
                                            pt_2d_list.push_back(
                                                _object_list_mm[m][tracer_id_match[m]].GetCenterPos()
                                            );
                                        } 
                                    }

                                    Matrix<double> pt_world(3,1);
                                    double error_3d=0.0;

                                    Triangulation(
                                        cam_cur_id+1,
                                        pt_2d_list,
                                        pt_world,
                                        error_3d
                                    );

                                    if (error_3d > _tor_3d)
                                    {
                                        break;
                                    }

                                    if (_check_cam_id < _n_cam)
                                    {
                                        int cam_next_id = cam_cur_id + 1;

                                        CheckTracerMatch (
                                            cam_next_id,
                                            tracer_id_match_list,
                                            tracer_id_match, 
                                            pt_world,
                                            error_list
                                        );

                                        tracer_id_match.pop_back();
                                    }
                                    else 
                                    {
                                        tracer_id_match_list.push_back(tracer_id_match);

                                        error_list.push_back(error_3d);

                                        tracer_id_match.pop_back();
                                    }
                                }
                            }
                            
                        }
                    }
                }
            }
        }
        
        #pragma endregion 
    }
    //                                  //
    // cam_cur_id = 2 ~ _check_cam_id-1 //
    //                                  //
    else
    {
        std::vector<Matrix<double>> line_of_sight_list;
        std::vector<Matrix<double>> pt_list;

        // Line of sight
        for (int i = 0; i < cam_cur_id; i ++)
        {       
            // project from cam_i onto cam_cur
            // find out tracer position in camera coordinate [mm]
            // find out tracer image position in world coordinate
            // _cam_list[i].ImgMMToWorld(_object_list_mm[i][tracer_id_match[i]].GetCenterPos())

            // position of cam_i center on cam_cur
            // Matrix<double> pt_center_i_on_cur = _cam_list[cam_cur_id].WorldToImgMM(_cam_list[i].GetCenterInWorld());
            pt_list.push_back(
                _cam_list[cam_cur_id].WorldToImgMM(_cam_list[i].GetCenterInWorld())
            );
            // position of tracer center on cam_i projected on cam_cur
            // Matrix<double> pt_tracer_i_on_cur = _cam_list[cam_cur_id].WorldToImgMM(pt_world);

            // unit vector in projected line of sight direction 
            Matrix<double> line_of_sight(
                _cam_list[cam_cur_id].WorldToImgMM(_cam_list[i].GetCenterInWorld()), 
                _cam_list[cam_cur_id].WorldToImgMM(
                    _cam_list[i].ImgMMToWorld(
                        _object_list_mm[i][tracer_id_match[i]].GetCenterPos()
                    )
                )
            );

            line_of_sight_list.push_back(line_of_sight);
        }
        PixelRange search_region_all = FindSearchRegion(cam_cur_id, pt_list, line_of_sight_list);
        // PixelRange search_region_all = FindSearchRegionOrig(cam_cur_id, pt_list, line_of_sight_list);

        int row_start = std::min(std::max(search_region_all._row_min, 0), n_row-1);
        int row_end   = std::min(std::max(search_region_all._row_max, 1), n_row);
        int col_start = std::min(std::max(search_region_all._col_min, 0), n_col-1);
        int col_end   = std::min(std::max(search_region_all._col_max, 1), n_col);

        // std::cout << "row:" << row_start << "," << row_end << ";"
        //           << "col:" << col_start << "," << col_end
        //           << "\n";

        for (int i = row_start; i < row_end; i ++)
        {
            for (int j = col_start; j < col_end; j ++)
            {
                // judge whether the distances between the candidate
                // and all the lines are all within the range 
                for (
                    int k = 0; 
                    k < _object_id_map[cam_cur_id][i][j].size(); 
                    k ++
                )
                {
                    int tracer_id;
                    tracer_id = _object_id_map[cam_cur_id][i][j][k];

                    if (tracer_id == -1)
                    {
                        break;
                    }

                    #pragma region check whether within error line 
                    
                    Matrix<double> candidate_2d (_object_list_mm[cam_cur_id][tracer_id].GetCenterPos());
                    bool in_range = true;
                    double dist = 0;
                    for (int m = 0; m < cam_cur_id; m ++)
                    {
                        dist = candidate_2d.Dist(pt_list[m], line_of_sight_list[m]);
                        if (dist > _tor_2d)
                        {
                            in_range = false;
                            break;
                        }
                    } 

                    // reproject onto the previous cam, then is within error line
                    if (in_range)
                    {
                        in_range = CheckReProject(cam_cur_id, tracer_id_match, candidate_2d);
                    }

                    #pragma endregion

                    if (in_range) 
                    {
                        tracer_id_match.push_back(tracer_id);

                        // if there is still other cameras to check
                        if (cam_cur_id >= 1 && cam_cur_id < _check_cam_id-1) 
                        {
                            int cam_next_id = cam_cur_id + 1;

                            // move onto the next camera and search candidates
                            FindTracerMatchOrig (
                                cam_next_id,
                                tracer_id_match_list,
                                tracer_id_match, 
                                error_list                  
                            );

                            tracer_id_match.pop_back(); // clear the match for the next candidate.
                        }
                        // if the current camera is the last one, then finalize the match.
                        else
                        {
                            //test 3d distance by triangulation                 
                            std::vector<Matrix<double>> pt_2d_list;
                            for (int m = 0; m < cam_cur_id+1; m ++)
                            {
                                if (m == cam_cur_id)
                                {
                                    pt_2d_list.push_back(
                                        _object_list_mm[m][tracer_id].GetCenterPos()
                                    ); 
                                }
                                else 
                                {
                                    pt_2d_list.push_back(
                                        _object_list_mm[m][tracer_id_match[m]].GetCenterPos()
                                    );
                                } 
                            }

                            Matrix<double> pt_world(3,1);
                            double error_3d=0.0;

                            Triangulation(
                                cam_cur_id+1,
                                pt_2d_list,
                                pt_world,
                                error_3d
                            );

                            if (error_3d > _tor_3d)
                            {
                                tracer_id_match.pop_back();
                            }
                            else
                            {
                                if (_check_cam_id < _n_cam)
                                {
                                    int cam_next_id = cam_cur_id + 1;

                                    CheckTracerMatch (
                                        cam_next_id,
                                        tracer_id_match_list,
                                        tracer_id_match, 
                                        pt_world,
                                        error_list
                                    );

                                    tracer_id_match.pop_back();
                                }
                                else 
                                {
                                    // // Print info 
                                    // std::cout << "Newly found match list: ";
                                    // for (int id = 0; id < _n_cam; id ++)
                                    // {
                                    //     std::cout << tracer_id_match[id] 
                                    //               << " ";
                                    // }
                                    // std::cout << std::endl;

                                    tracer_id_match_list.push_back(tracer_id_match);

                                    error_list.push_back(error_3d);

                                    tracer_id_match.pop_back();
                                }
                            }
                            
                        }
                    }
                }  
            }
        }
    }
}


template<class T>
bool StereoMatch<T>::CheckReProject (
    int cam_cur_id, 
    std::vector<int> const& tracer_id_match, 
    Matrix<double> const& pt_mm
)
{
    double dist = 0;
    Matrix<double> temp(3,1);
    for (int i = 0; i < cam_cur_id; i ++)
    {
        temp = _cam_list[i].WorldToImgMM(
            _cam_list[cam_cur_id].ImgMMToWorld(pt_mm)
        );
        Matrix<double> unit(
            _cam_list[i].WorldToImgMM(_cam_list[cam_cur_id].GetCenterInWorld()),
            temp
        );
        dist = _object_list_mm[i][tracer_id_match[i]].GetCenterPos().Dist(temp, unit);

        if (dist > _tor_2d)
        {
            return false;
        }
    }
    return true;
}


template<class T>
void StereoMatch<T>::CheckTracerMatch (
    int cam_cur_id, 
    std::deque<std::vector<int>>& tracer_id_match_list,
    std::vector<int>& tracer_id_match, 
    Matrix<double>& pt_world,
    std::deque<double>& error_list
)
{
    int n_row = _cam_list[cam_cur_id].GetNpixh();
    int n_col = _cam_list[cam_cur_id].GetNpixw();

    if (cam_cur_id < 1)
    {
        std::cerr << "StereoMatch error: "
                  << "camera id is too small to "
                  << "start checking tracer match!"
                  << std::endl;
        throw;
    }


    std::vector<Matrix<double>> line_of_sight_list;
    std::vector<Matrix<double>> pt_list;
    // Line of sight
    for (int i = 0; i < cam_cur_id; i ++)
    {       
        // project from cam_i onto cam_cur
        // find out tracer position in camera coordinate [mm]
        // find out tracer image position in world coordinate
        // _cam_list[i].ImgMMToWorld(_object_list_mm[i][tracer_id_match[i]].GetCenterPos())

        // position of cam_i center on cam_cur
        // Matrix<double> pt_center_i_on_cur = _cam_list[cam_cur_id].WorldToImgMM(_cam_list[i].GetCenterInWorld());
        pt_list.push_back(
            _cam_list[cam_cur_id].WorldToImgMM(_cam_list[i].GetCenterInWorld())
        );
        // position of tracer center on cam_i projected on cam_cur
        // Matrix<double> pt_tracer_i_on_cur = _cam_list[cam_cur_id].WorldToImgMM(pt_world);

        // unit vector in projected line of sight direction 
        Matrix<double> line_of_sight(
            _cam_list[cam_cur_id].WorldToImgMM(_cam_list[i].GetCenterInWorld()), 
            _cam_list[cam_cur_id].WorldToImgMM(
                _cam_list[i].ImgMMToWorld(
                    _object_list_mm[i][tracer_id_match[i]].GetCenterPos()
                )
            )
        );

        line_of_sight_list.push_back(line_of_sight);
    }


    // Project pt_world onto testing cam[cam_cur_id]
    double factor = 1.0;
    Matrix<double> pt_pixel_upright (
        _cam_list[cam_cur_id].ImgMMToPixel(
            _cam_list[cam_cur_id].WorldToImgMM(pt_world)
            + Matrix<double>(3,1, factor*_tor_2d, factor*_tor_2d, 0.0)
        )
    );
    Matrix<double> pt_pixel_lowleft (
        _cam_list[cam_cur_id].ImgMMToPixel(
            _cam_list[cam_cur_id].WorldToImgMM(pt_world)
            - Matrix<double>(3,1, factor*_tor_2d, factor*_tor_2d, 0.0)
        )
    );
    //Matrix<double> pt_pixel 
    //= _cam_list[cam_cur_id].ImgMMToPixel(
    //    _cam_list[cam_cur_id].WorldToImgMM(pt_world)
    //);

    int row_start = int( 
        std::min(
            std::max(
                int(pt_pixel_upright(1,0))
                , 0
            ), 
            n_row - 1 
        ) 
    );
    int row_end   = int( 
        std::min(
            std::max(
                int(pt_pixel_lowleft(1,0)+1.0)
                , 1
            ), 
            n_row
        ) 
    );
    int col_start = int( 
        std::min(
            std::max(
                int(pt_pixel_lowleft(0,0))
                , 0
            ), 
            n_col - 1 
        ) 
    );
    int col_end   = int( 
        std::min(
            std::max(
                int(pt_pixel_upright(0,0)+1.0)
                , 1
            ), 
            n_col
        ) 
    );

    // Search around pt_pixel for candidates
    for (int i = row_start; i < row_end; i ++)
    {
        for (int j = col_start; j < col_end; j ++)
        {
            for (
                int k = 0; 
                k < _object_id_map[cam_cur_id][i][j].size(); 
                k ++
            )
            {
                int tracer_id;
                tracer_id = _object_id_map[cam_cur_id][i][j][k];

                if (tracer_id == -1)
                {
                    break;
                }

                #pragma region check whether within error line 
                    
                Matrix<double> candidate_2d (_object_list_mm[cam_cur_id][tracer_id].GetCenterPos());
                bool in_range = true;
                double dist = 0;
                for (int m = 0; m < cam_cur_id; m ++)
                {
                    dist = candidate_2d.Dist(pt_list[m], line_of_sight_list[m]);
                    if (dist > _tor_2d)
                    {
                        in_range = false;
                        break;
                    }
                } 

                // reproject onto the previous cam, then is within error line
                if (in_range)
                {
                    in_range = CheckReProject(cam_cur_id, tracer_id_match, candidate_2d);
                }

                #pragma endregion

                // check 3D here, triangulation
                if (in_range) 
                {
                    tracer_id_match.push_back(tracer_id);

                    // if there is still other cameras to check
                    if (cam_cur_id >= 1 && cam_cur_id < _n_cam-1) 
                    {
                        int cam_next_id = cam_cur_id + 1;

                        // move onto the next camera and search candidates
                        CheckTracerMatch (
                            cam_next_id, 
                            tracer_id_match_list,
                            tracer_id_match, 
                            pt_world,
                            error_list
                        );
                        tracer_id_match.pop_back(); // clear the match for the next candidate.
                    }
                    // if the current camera is the last one, then finalize the match.
                    else
                    {
                        // // Print info 
                        // std::cout << "Newly found match list: ";
                        // for (int id = 0; id < _n_cam; id ++)
                        // {
                        //     std::cout << tracer_id_match[id] 
                        //               << " ";
                        // }
                        // std::cout << "\n";

                        //test 3d distance by triangulation                 
                        std::vector<Matrix<double>> pt_2d_list;
                        for (int m = 0; m < cam_cur_id+1; m ++)
                        {
                            if (m == cam_cur_id)
                            {
                                pt_2d_list.push_back(
                                    _object_list_mm[m][tracer_id].GetCenterPos()
                                ); 
                            }
                            else 
                            {
                                pt_2d_list.push_back(
                                    _object_list_mm[m][tracer_id_match[m]].GetCenterPos()
                                );
                            } 
                        }

                        Matrix<double> pt_world(3,1);
                        double error_3d=0.0;

                        Triangulation(
                            cam_cur_id+1,
                            pt_2d_list,
                            pt_world,
                            error_3d
                        );

                        
                        if (error_3d <= _tor_3d)
                        {
                            tracer_id_match_list.push_back(tracer_id_match);

                            error_list.push_back(error_3d);
                        }   

                        tracer_id_match.pop_back();
                    }
                }
            }
        }
    }
}


template<class T>
void StereoMatch<T>::DeleteGohstTracerMatch (std::vector<std::vector<TracerInfo>>& object_list_pixel)
{
    std::cout << "Start deleting gohst match!" << std::endl;

    std::vector<std::vector<int>> object_id_match_map; // indicate whether the 2D point is used
    // object_id_match_map.reserve(2 + _n_cam);

    int temp = 0;
    while (temp < 5)
    {
        try 
        {
            object_id_match_map.reserve(2 + _n_cam);

            break;
        }
        catch (...)
        {    
            if (temp == 4)
            {
                std::runtime_error("DeleteGohstTracerMatch: object_id_match_map no enough space");
                throw;
            }
        }
        temp += 1;
    }

    // initialization
    for (int cam_id = 0; cam_id < _n_cam; cam_id ++)
    {
        std::vector<int> object_id_match_map_one;

        int temp = 0;
        while (temp < 5)
        {
            try 
            {
                object_id_match_map_one.reserve(5 + object_list_pixel[cam_id].size());

                break;
            }
            catch (...)
            {    
                if (temp == 4)
                {
                    std::runtime_error("DeleteGohstTracerMatch: object_id_match_map_one no enough space");
                    throw;
                }
            }
            temp += 1;
        }
        // object_id_match_map_one.reserve(5 + object_list_pixel[cam_id].size());

        object_id_match_map_one.resize(object_list_pixel[cam_id].size(), 0);
        object_id_match_map.push_back(object_id_match_map_one);
    }

    int n = _error_list.size();
    std::vector<int> ranked_index(n);
    myMATH::MergeSort<double>(ranked_index, _error_list);
    for (int i = 0; i < n; i ++)
    {
        int sum = 0; // check whether all 2d pts in a match is available.
        for (int cam_id = 0; cam_id < _n_cam; cam_id ++)
        {
            sum += object_id_match_map[cam_id][_object_id_match_list[ranked_index[i]][cam_id]];
        }
        if (! sum) // if available
        {
            std::vector<Matrix<double>> pt_2d_list;

            int temp = 0;
            while (temp < 5)
            {
                try 
                {
                    pt_2d_list.reserve(_n_cam + 2);

                    break;
                }
                catch (...)
                {    
                    if (temp == 4)
                    {
                        std::runtime_error("DeleteGohstTracerMatch: pt_2d_list no enough space");
                        throw;
                    }
                }
                temp += 1;
            }

            double error;
            Matrix<double> pt_world(3,1); // 3D pos
            TracerInfo tracer;
            for (int cam_id = 0; cam_id < _n_cam; cam_id ++)
            {
                int tracer_id = _object_id_match_list[ranked_index[i]][cam_id];
                object_id_match_map[cam_id][tracer_id] = 1; // indicate the 2D pt is used
                pt_2d_list.push_back(_object_list_mm[cam_id][tracer_id].GetCenterPos());
                tracer.AddMatchInfo(object_list_pixel[cam_id][tracer_id]);
            }
            Triangulation(_n_cam, pt_2d_list, pt_world, error);
            tracer.SetError(error);
            tracer.SetCenterPos(pt_world);

            _object_info_match_list.push_back(tracer);
        }
    }
    
    _n_after_del = _object_info_match_list.size();
    _n_del = _n_before_del - _n_after_del;
    std::cout << "Finish deleting gohst match!" << std::endl;
    std::cout << "n_del = " << _n_del << ", "
              << "n_after_del = " << _n_after_del << "."
              << std::endl;
}


template<class T>
void StereoMatch<T>::DeleteGohstTracerMatchNew (std::vector<std::vector<TracerInfo>>& object_list_pixel)
{
    std::cout << "Start deleting gohst match!" << std::endl;

    int n_match = _object_id_match_list.size();

    std::vector<double> object_id_match_rank;
    std::vector<std::vector<int>> object_id_match_num;
    std::vector<std::vector<int>> object_id_match_map; // indicate whether the 2D point is used
    // object_id_match_map.reserve(2 + _n_cam);

    int temp = 0;
    while (temp < 5)
    {
        try 
        {
            object_id_match_map.reserve(2 + _n_cam);

            break;
        }
        catch (...)
        {    
            if (temp == 4)
            {
                std::runtime_error("DeleteGohstTracerMatch: object_id_match_map no enough space");
                throw;
            }
        }
        temp += 1;
    }
    
    std::cout << "1" << std::endl;

    // initialization
    for (int cam_id = 0; cam_id < _n_cam; cam_id ++)
    {
        std::vector<int> object_id_match_map_one;

        int temp = 0;
        while (temp < 5)
        {
            try 
            {
                object_id_match_map_one.reserve(5 + object_list_pixel[cam_id].size());

                break;
            }
            catch (...)
            {    
                if (temp == 4)
                {
                    std::runtime_error("DeleteGohstTracerMatch: object_id_match_map_one no enough space");
                    throw;
                }
            }
            temp += 1;
        }
        // object_id_match_map_one.reserve(5 + object_list_pixel[cam_id].size());

        object_id_match_map_one.resize(object_list_pixel[cam_id].size(), 0);
        
        object_id_match_num.push_back(object_id_match_map_one);
        object_id_match_map.push_back(object_id_match_map_one);
    }

    std::cout << "2" << std::endl;

    // calculate number of each id used (might be parallel)
    for (int i = 0; i < n_match; i ++)
    {
        if (_object_id_match_list[i].size() != _n_cam)
        {
            std::cout << "wrong: " << _object_id_match_list[i].size() << std::endl;
            // throw;
            continue;
        }
        for (int j = 0; j < _n_cam; j ++)
        {
            int id = _object_id_match_list[i][j];
            // object_id_match_num[j][id] += 1;
            object_id_match_num.at(j).at(id) += 1;
        }  
    }

    std::cout << "3" << std::endl;

    // calculate rank of each match pair (might be parallel)
    // if (_n_thread != 0)
    // {
    //     omp_set_num_threads(_n_thread);
    // }
    // #pragma omp parallel
    {
        // #pragma omp for
        for (int i = 0; i < n_match; i ++)
        {
            // #pragma critical
            {
                double value = 5.175e2 * _error_list[i];
                // double value = 5.175e2 / _error_list[i];
                for (int j = 0; j < _n_cam; j ++)
                {
                    int id = _object_id_match_list[i][j];
                    value += double(object_id_match_num[j][id]);
                }
            
                object_id_match_rank.push_back(value);
            }
            
        }
    }

    std::cout << "4" << std::endl;

    std::vector<int> ranked_index (n_match);
    myMATH::MergeSort<double>(ranked_index, object_id_match_rank);

    std::cout << "5" << std::endl;

    for (int i = 0; i < n_match; i ++)
    // for (int i = n_match-1; i >-1; i --)
    {
        int sum = 0; // check whether all 2d pts in a match is available.
        for (int cam_id = 0; cam_id < _n_cam; cam_id ++)
        {
            sum += object_id_match_map[cam_id][_object_id_match_list[ranked_index[i]][cam_id]];
        }
        if (! sum) // if available
        {
            std::vector<Matrix<double>> pt_2d_list;

            int temp = 0;
            while (temp < 5)
            {
                try 
                {
                    pt_2d_list.reserve(_n_cam + 2);

                    break;
                }
                catch (...)
                {    
                    if (temp == 4)
                    {
                        std::runtime_error("DeleteGohstTracerMatch: pt_2d_list no enough space");
                        throw;
                    }
                }
                temp += 1;
            }

            double error;
            Matrix<double> pt_world(3,1); // 3D pos
            TracerInfo tracer;
            for (int cam_id = 0; cam_id < _n_cam; cam_id ++)
            {
                int tracer_id = _object_id_match_list[ranked_index[i]][cam_id];
                object_id_match_map[cam_id][tracer_id] = 1; // indicate the 2D pt is used
                pt_2d_list.push_back(_object_list_mm[cam_id][tracer_id].GetCenterPos());
                tracer.AddMatchInfo(object_list_pixel[cam_id][tracer_id]);
            }
            Triangulation(_n_cam, pt_2d_list, pt_world, error);
            tracer.SetError(error);
            tracer.SetCenterPos(pt_world);

            _object_info_match_list.push_back(tracer);
        }
    }
    
    _n_after_del = _object_info_match_list.size();
    _n_del = _n_before_del - _n_after_del;
    std::cout << "Finish deleting gohst match!" << std::endl;
    std::cout << "n_del = " << _n_del << ", "
              << "n_after_del = " << _n_after_del << "."
              << std::endl;
}


template<class T>
void StereoMatch<T>::Merge (std::vector<int>& A, int id_begin, int id_middle, int id_end, std::vector<int>& B, std::vector<int>& C)
{
    int i = id_begin;
    int j = id_middle;

    for (int k = id_begin; k < id_end; k ++)
    {
        if (i < id_middle && j >= id_end)
        {
            B[k] = A[i];
            i ++;
        }
        else if (i < id_middle && C[A[i]] > C[A[j]])
        {
            B[k] = A[i];
            i ++;
        }
        else if (i < id_middle && C[A[i]] == C[A[j]])
        {
            // when the preference numbers are the same
            // then compare the error
            if (_error_list[A[i]] < _error_list[A[j]])
            {
                B[k] = A[i];
                i ++;
            }
            else 
            {
                B[k] = A[j];
                j ++; 
            }
        }
        else
        {
            B[k] = A[j];
            j ++; 
        }
    }
}


template<class T>
void StereoMatch<T>::SplitMerge (std::vector<int>& B, int id_begin, int id_end, std::vector<int>& A, std::vector<int>& C)
{
    if (id_end - id_begin < 2)
    {
        return;
    }
    int id_middle = (id_begin + id_end) / 2;
    SplitMerge(A, id_begin, id_middle, B, C);
    SplitMerge(A, id_middle, id_end, B, C);
    Merge(B, id_begin, id_middle, id_end, A, C);
}


template<class T>
void StereoMatch<T>::SortPreference (std::vector<int>& sorted_id_list, std::vector<int>& object_id_match_rank)
{
    int n_match = _object_id_match_list.size();

    std::vector<int> B;
    B.reserve(n_match);
    for (int i = 0; i < n_match; i ++)
    {
        B[i] = i;
    }

    SplitMerge(B, 0, n_match, sorted_id_list, object_id_match_rank);
}


template<class T>
void StereoMatch<T>::DeleteGohstTracerMatchOrig (std::vector<std::vector<TracerInfo>>& object_list_pixel)
{
    std::cout << "Start deleting gohst match!" << std::endl;

    int n_match = _object_id_match_list.size();

    //------------------------------------//
    // update object_id_match_rank
    std::vector<int> object_id_match_rank;
    object_id_match_rank.reserve(n_match);
    object_id_match_rank.resize(n_match, 0);
    for (int i = 0; i < _n_cam; i ++)
    {
        int n_2D_object = object_list_pixel[i].size();

        // buffer_id_map: save the match id with smallest error 
        std::vector<int> buffer_id_map;
        buffer_id_map.reserve(n_2D_object);
        buffer_id_map.resize(n_2D_object, -1);
        for (int j = 0; j < n_match; j ++)
        {
            int match_id_2D = _object_id_match_list[j][i];
            
            if (buffer_id_map[match_id_2D] == -1) 
            {
                buffer_id_map[match_id_2D] = j;
            }
            else // there is a object_id filling
            {
                int object_id = buffer_id_map[match_id_2D];
                if (_error_list[object_id] > _error_list[j])
                {
                    // replaced by a object_id with smaller error 
                    buffer_id_map[match_id_2D] = j;
                }
            }
        }

        // update the rank
        for (int j = 0; j < n_2D_object; j ++)
        {
            if (buffer_id_map[j] == -1)
            {
                continue;
            }
            int object_id = buffer_id_map[j];
            object_id_match_rank[object_id] += 1;
        }
    }

    // sort object_id_match_rank based on rank 
    // merge sort
    // small to large 
    std::vector<int> sorted_id_list;
    sorted_id_list.reserve(n_match);
    for (int i = 0; i < n_match; i ++)
    {
        sorted_id_list[i] = i;
    }
    SortPreference(sorted_id_list, object_id_match_rank);

    // map 2D particle id to match id with smallest error
    std::vector<std::vector<int>> object_id_match_map;
    for (int i = 0; i < _n_cam; i ++)
    {
        int n_2D_object = object_list_pixel[i].size();

        std::vector<int> object_id_match_map_one;
        object_id_match_map_one.reserve(n_2D_object);
        object_id_match_map_one.resize(n_2D_object, -1);

        object_id_match_map.push_back(object_id_match_map_one);
    }

    std::vector<bool> is_select;
    is_select.reserve(n_match);
    is_select.resize(n_match, false); 

    for (int i = 0; i < n_match; i ++)
    {
        // sorted object_id (3D) based on rank and error
        int object_id = sorted_id_list[i];  

        int n_equal = 0;
        double sum_error = 0.0;

        int is_rank_small = 0;
        for (int j = 0; j < _n_cam; j ++)
        {
            int match_id_2D = _object_id_match_list[object_id][j];
            int map_id = object_id_match_map[j][match_id_2D];

            if (map_id == -1)
            {
                // when there is no match for the checked particle
                continue;          
            }
            else if (object_id_match_rank[object_id] < object_id_match_rank[map_id] && is_rank_small==0) 
            {
                is_rank_small = 1;
            }
            else if (object_id_match_rank[object_id] == object_id_match_rank[map_id])
            {
                // when the preference # is the same as the current preference # of those filled matches,
                // then we should calculate the average of the preference error of those filled matches
                n_equal ++;
                sum_error += _error_list[map_id];
            }
        }
        if (
            ( 
              n_equal > 0 && 
              _error_list[object_id] < sum_error / n_equal
            ) || n_equal == 0 && is_rank_small == 0
        )
        {
            // 1) if there is same preference number (for certain cam),
            // and the error is less than the averaged of those filled match, 
            // then replace all the old match with the new one
            // 2) no repeating

            for (int j = 0; j < _n_cam; j ++)
            {
                int match_id_2D = _object_id_match_list[object_id][j];
                int map_id = object_id_match_map[j][match_id_2D];

                if (map_id != -1)
                {
                    is_select[map_id] = false;
                    for (int k = 0; k < _n_cam; k ++)
                    {
                        int origin_match_id_2D = _object_id_match_list[map_id][k];

                        object_id_match_map[k][origin_match_id_2D] = -1;
                    }
                }
                object_id_match_map[j][match_id_2D] = object_id;
            }

            is_select[object_id] = true;
        }
    }
    
    // collect new match list 
    for (int i = 0; i < n_match; i ++)
    {
        
        if (is_select[i])
        {
            std::vector<Matrix<double>> pt_2d_list;
            pt_2d_list.reserve(_n_cam);
            
            double error;
            Matrix<double> pt_world(3,1); // 3D pos
            TracerInfo tracer;
            for (int cam_id = 0; cam_id < _n_cam; cam_id ++)
            {
                int tracer_id = _object_id_match_list[i][cam_id];

                pt_2d_list.push_back(_object_list_mm[cam_id][tracer_id].GetCenterPos());
                
                tracer.AddMatchInfo(object_list_pixel[cam_id][tracer_id]);
            }
            Triangulation(_n_cam, pt_2d_list, pt_world, error);
            tracer.SetError(error);
            tracer.SetCenterPos(pt_world);

            _object_info_match_list.push_back(tracer);
        }
    }

    //------------------------------------//
    _n_after_del = _object_info_match_list.size();
    _n_del = _n_before_del - _n_after_del;
    std::cout << "Finish deleting gohst match!" << std::endl;
    std::cout << "n_del = " << _n_del << ", "
              << "n_after_del = " << _n_after_del << "."
              << std::endl;
}




template<class T>
void StereoMatch<T>::FillObjectInfo (std::vector<std::vector<T>>& object_list_pixel)
{
    for (int i = 0; i < _object_id_match_list.size(); i ++)
    {
        T object;
        std::vector<Matrix<double>> pt_2d_list;
        double error;
        Matrix<double> pt_world(3,1);

        for (int cam_id = 0; cam_id < _n_cam; cam_id ++)
        {
            int object_id = _object_id_match_list[i][cam_id];
            pt_2d_list.push_back(_object_list_mm[cam_id][object_id].GetCenterPos());
            object.AddMatchInfo(object_list_pixel[cam_id][object_id]);
        }
        Triangulation(_n_cam, pt_2d_list, pt_world, error);
        object.SetError(error);
        object.SetCenterPos(pt_world);

        _object_info_match_list.push_back(object);
    }
}


// Return list of tracer_id
//
//
template<class T>
void StereoMatch<T>::TracerMatch (std::vector<std::vector<TracerInfo>>& object_list_pixel)
{
    std::cout << "Tracer match start!" << std::endl;

    if (_n_cam < 2)
    {
        std::cerr << "StereoMatch error: "
                  << "Need at least 2 cameras for matching" 
                  << std::endl;
        throw error_size;
    }
    if (int(object_list_pixel.size()) != _n_cam)
    {
        std::cerr << "The size of tracer list is: " 
                  << object_list_pixel.size() 
                  << ", is different from the size of camera list: "
                  << _n_cam << "."
                  << std::endl;
        throw error_size;
    }

    SetObjectListMMFromPixel (object_list_pixel);
    MakeObjectIDMap (object_list_pixel);

    // for 1st camera, draw a line of sight through each particle on its image plane.
    // project these lines of sight onto the image planes of 2nd camera.
    // particles within mindist_2D of these lines are candidate matches between first 2 cams.
    // then project 2 line of sights from each particle pair of 1st n 2nd cam onto 3rd cam.
    // particles within mindist_1D of the intersection point are candidate matches from 3rd cam.
    // repeat similarly for subsequent cams

    int num_object = object_list_pixel[0].size();

    if (_n_thread != 0)
    {
        omp_set_num_threads(_n_thread);
    }
    #pragma omp parallel
    {
        #pragma omp for 
        for (int tracer_id = 0; tracer_id < _object_list_mm[0].size(); tracer_id ++)
        {
            std::deque<std::vector<int>> tracer_id_match_list; // match list for the particle i in the first camera

            std::vector<int> tracer_id_match; 
            tracer_id_match.reserve(_n_cam);
            tracer_id_match.push_back(tracer_id);

            std::deque<double> error_list;
            
            // zsj debug, 2022/12/06
            // FindTracerMatch (
            //     1,
            //     tracer_id_match_list,
            //     tracer_id_match,
            //     error_list
            // );
            FindTracerMatchOrig (
                1,
                tracer_id_match_list,
                tracer_id_match,
                error_list
            );

            #pragma omp critical
            {   
                if (error_list.size() != tracer_id_match_list.size())
                {
                    std::cout << "error_list.size = " << error_list.size()
                              << "tracer_id_match_list.size = " << tracer_id_match_list.size() << "\n";
                    throw;
                }
                for (int i = 0; i < tracer_id_match_list.size(); i ++) 
                {      
                    if (tracer_id_match_list[i].size() != _n_cam)
                    {
                        continue;
                    }

                    int n_error_size_pre = _error_list.size();
                    int n_match_size_pre = _object_id_match_list.size();

                    _error_list.push_back(error_list[i]);
                    _object_id_match_list.push_back(tracer_id_match_list[i]);  
                    
                    int n_error_size = _error_list.size();
                    int n_match_size = _object_id_match_list.size();

                    if (n_error_size == n_error_size_pre || n_match_size == n_match_size_pre)
                    {
                        std::cout << "Go into trouble!" << "\n"
                                  << "_error_list.size != _object_id_match_list.size; "
                                  << "_error_list.size = " << n_error_size
                                << "," << n_error_size_pre << ", _object_id_match_list.size = " << n_match_size << "," << n_match_size_pre
                                << "\n";

                        throw;
                    }
                  
                }
            }
        }
    }
    
    if (_error_list.size() != _object_id_match_list.size())
    {
        std::cout << "Go into trouble!" << "\n"
                  << "_error_list.size != _object_id_match_list.size; "
                  << "_error_list.size = " << _error_list.size()
                  << ", _object_id_match_list.size = " << _object_id_match_list.size()
                  << "\n";
        throw;
    }

    _n_before_del = _object_id_match_list.size();
    std::cout << "Finish stereomatching!" << std::endl;
    std::cout << "n_before_del = " << _n_before_del << "."
              << std::endl;
}


template<class T>
void StereoMatch<T>::Match (
    std::vector<std::vector<T>>& object_list_pixel, 
    int is_delete_ghost
)
{
    if (typeid(T) == typeid(TracerInfo))
    {   
        TracerMatch(object_list_pixel);
        if (is_delete_ghost)
        {   
            // debug zsj 01/04/2023
            // DeleteGohstTracerMatch(object_list_pixel);
            // DeleteGohstTracerMatchNew(object_list_pixel);
            DeleteGohstTracerMatchOrig(object_list_pixel);
        }
        else
        {
            FillObjectInfo(object_list_pixel);
        }
    }
    else 
    {
        std::cerr << "StereoMatch: " 
                  << "class " << typeid(T).name()
                  << "is not included in StereoMatch!"
                  << std::endl;
        throw;
    }
}


#endif