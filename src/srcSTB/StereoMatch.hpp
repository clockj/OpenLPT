#ifndef STEREOMATCH_HPP
#define STEREOMATCH_HPP

#include "StereoMatch.h"

StereoMatch::StereoMatch(StereoMatchParam param, CamList& cam_list) 
    : _param(param), _cam_list(cam_list), _n_cam_use(cam_list.useid_list.size()) {}

template<class T>
void StereoMatch::match(std::vector<T>& obj3d_list, std::vector<std::vector<T>> const& obj2d_list)
{
    if (typeid(T) == typeid(Tracer2D))
    {   
        if (obj3d_list.size() > 0)
        {
            obj3d_list.clear();
        }

        tracerMatch(obj3d_list, obj2d_list);

        if (_param.is_delete_ghost)
        {   
            // debug zsj 01/04/2023
            // DeleteGohstTracerMatch(object_list_pixel);
            // DeleteGohstTracerMatchNew(object_list_pixel);
            deleteGohstTracerMatchOrig(obj2d_list);
        }
        else
        {
            fillObjectInfo(obj2d_list);
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

template<class T>
void StereoMatch::createObjIDMap (std::vector<std::vector<T>> const& obj2d_list)
{
    if (_objID_map_list.size() > 0)
    {
        _objID_map_list.clear();
    }

    int row_id, col_id;
    int cam_id;
    for (int i = 0; i < _n_cam_use; i ++)
    {
        cam_id = _cam_list.useid_list[i];

        ObjIDMap objID_map;
        objID_map.config(
            _cam_list.cam_list[cam_id].getNRow(), 
            _cam_list.cam_list[cam_id].getNCol()
        );
        
        for (int j = 0; j < obj2d_list[i].size(); j ++)
        {
            row_id = obj2d_list[i][j]._pt_center[1]; // img_y
            col_id = obj2d_list[i][j]._pt_center[0]; // img_x

            if (objID_map(row_id, col_id)[0] == -1)
            {
                objID_map(row_id, col_id)[0] = j;
            }
            else
            {
                objID_map(row_id, col_id).push_back(j);
            }
        }

        _objID_map_list.push_back(objID_map);
    }
}


//              //
// Tracer match //
//              //

// Return list of tracer_id
void StereoMatch::tracerMatch (std::vector<Tracer3D>& tr3d_list, std::vector<std::vector<Tracer2D>> const& tr2d_list)
{
    // std::cout << "\t Tracer match start!" << std::endl;

    if (_n_cam_use < 2)
    {
        std::cerr << "StereoMatch::tracerMatch error at line " << __LINE__ << ": \n"
                  << "Need at least 2 cameras for matching" 
                  << std::endl;
        throw error_size;
    }
    if (int(tr2d_list.size()) != _n_cam_use)
    {
        std::cerr << "StereoMatch::tracerMatch error at line " << __LINE__ << ":\n"
                  << "The size of tracer list is: " 
                  << tr2d_list.size() 
                  << ", is different from the size of camera list: "
                  << _n_cam_use << "."
                  << std::endl;
        throw error_size;
    }

    createObjIDMap (tr2d_list);

    if (_param.n_thread != 0)
    {
        omp_set_num_threads(_param.n_thread);
    }
    #pragma omp parallel
    {
        // for 1st used camera, draw a line of sight through each particle on its image plane.
        // project these lines of sight onto the image planes of 2nd camera.
        // particles within mindist_2D of these lines are candidate matches between first 2 cams.
        // then project 2 line of sights from each particle pair of 1st n 2nd cam onto 3rd cam.
        // particles within mindist_1D of the intersection point are candidate matches from 3rd cam.
        // repeat similarly for subsequent cams
        #pragma omp for 
        for (int tr_id = 0; tr_id < tr2d_list[0].size(); tr_id ++)
        {
            std::deque<std::vector<int>> trID_match_list; // match list for the particle i in the first camera

            std::vector<int> trID_match; 
            trID_match.reserve(_n_cam_use);
            trID_match.push_back(tr_id);

            std::deque<double> error_list;
            
            findTracerMatch (
                1,                
                trID_match,
                trID_match_list,
                error_list,
                tr2d_list
            );

            #pragma omp critical
            {   
                if (error_list.size() != trID_match_list.size())
                {
                    std::cerr << "StereoMatch::tracerMatch error at line " << __LINE__ << ": \n" 
                              << "error_list.size = " << error_list.size()
                              << "tracer_id_match_list.size = " << trID_match_list.size() << "\n";
                    throw error_size;
                }
                for (int i = 0; i < trID_match_list.size(); i ++) 
                {      
                    if (trID_match_list[i].size() != _n_cam_use)
                    {
                        continue;
                    }

                    int n_error_size_pre = _error_list.size();
                    int n_match_size_pre = _objID_match_list.size();

                    _error_list.push_back(error_list[i]);
                    _objID_match_list.push_back(trID_match_list[i]);  
                    
                    int n_error_size = _error_list.size();
                    int n_match_size = _objID_match_list.size();
                    if (n_error_size == n_error_size_pre || n_match_size == n_match_size_pre)
                    {
                        std::cerr << "StereoMatch::tracerMatch error at line " << __LINE__ << ":\n"
                                  << "_error_list.size != _object_id_match_list.size; "
                                  << "_error_list.size = " << n_error_size
                                  << "," << n_error_size_pre << ", _object_id_match_list.size = " << n_match_size << "," << n_match_size_pre
                                  << "\n";
                        throw error_size;
                    }
                }
            }
        }
    }
    
    if (_error_list.size() != _objID_match_list.size())
    {
        std::cerr << "StereoMatch::tracerMatch error at line " << __LINE__ << ":\n"
                  << "_error_list.size != _objID_match_list.size; "
                  << "_error_list.size = " << _error_list.size()
                  << ", _object_id_match_list.size = " << _objID_match_list.size()
                  << "\n";
        throw error_size;
    }

    _n_before_del = _objID_match_list.size();
    std::cout << "\tFinish stereomatching: " 
              << "n_before_del = " << _n_before_del << "."
              << std::endl;
}


void StereoMatch::findTracerMatch (
    int id, // id = 1 => 2nd used cam 
    std::vector<int> trID_match,
    std::deque<std::vector<int>>& trID_match_list,
    std::deque<double>& error_list,
    std::vector<std::vector<Tracer2D>> const& tr2d_list
)
{
    int camID_curr = _cam_list.useid_list[id];
    int n_row = _cam_list.cam_list[camID_curr].getNRow();
    int n_col = _cam_list.cam_list[camID_curr].getNCol();

    // Line of sight
    std::vector<Line2D> sight2D_list;
    std::vector<Line3D> sight3D_list;
    Line2D sight2D;
    Line3D sight3D;
    Pt2D pt2d_1;
    Pt2D pt2d_2;
    Pt2D unit2d;
    for (int i = 0; i < id; i ++)
    {       
        // project from cam_prev onto cam_curr
        int camID_prev = _cam_list.useid_list[i];

        // get 3d light of sight from cam_prev
        sight3D = _cam_list.cam_list[camID_prev].lineOfSight(tr2d_list[i][trID_match[i]]._pt_center);
        sight3D_list.push_back(sight3D);

        // project 3d light of sight onto cam_curr (3d line => 2d line)
        pt2d_1 = _cam_list.cam_list[camID_curr].project(sight3D.pt);
        pt2d_2 = _cam_list.cam_list[camID_curr].project(sight3D.pt + sight3D.unit_vector);
        unit2d = myMATH::createUnitVector(pt2d_1, pt2d_2);
        sight2D.pt = pt2d_1;
        sight2D.unit_vector = unit2d;
        sight2D_list.push_back(sight2D);
    }
    sight3D_list.push_back(sight3D);


    //                       //
    // id = 1 (2nd used cam) //
    //                       //
    // tracer marker around the line of sight 
    // intialize a 2D map with value 0
    // optimization for cam_cur_id == 1
    Pt3D pt3d;
    if (id == 1)  
    {
        // find tracer candidates around the line of sight         
        // ObjectMarkerAroundLine   

        // calculate the parallel border lines: represented by the ref points
        Line2D sight2D_plus;
        Line2D sight2D_minus;
        // calculate the perpendicular unit vecotr
        unit2d[0] = sight2D.unit_vector[1];
        unit2d[1] = -sight2D.unit_vector[0];
        // calculate the ref points
        pt2d_1 = sight2D.pt + unit2d * _param.tor_2d;
        pt2d_2 = sight2D.pt - unit2d * _param.tor_2d;
        sight2D_plus.pt = pt2d_1;
        sight2D_plus.unit_vector = sight2D.unit_vector;
        sight2D_minus.pt = pt2d_2;
        sight2D_minus.unit_vector = sight2D.unit_vector;

        // Search the candidate within the torlerance near a line
        // Determine the border of the torlerance by iterating each x pixel or y pixel depending on the slope.
        // x_pixel (col id), y_pixel (row id)
        //  if the |slope| > 1, iterate every y_pixel (row)
        //  else, iterate every x_pixel (col) 
        Line2D sight2D_axis;
        if (std::fabs(sight2D.unit_vector[1]) > std::fabs(sight2D.unit_vector[0]))
        {
            sight2D_axis.unit_vector = Pt2D(1,0);
            int x_pixel_1, x_pixel_2, min, max;   

            int y_start = 0;
            int y_end = n_row;
            for (int y_pixel = y_start; y_pixel < y_end; y_pixel ++)
            {
                sight2D_axis.pt[0] = 0;
                sight2D_axis.pt[1] = y_pixel;

                // Get x_pixel (col id) from two lines 
                pt2d_1 = myMATH::crossPoint(sight2D_axis, sight2D_plus);
                pt2d_2 = myMATH::crossPoint(sight2D_axis, sight2D_minus);
                x_pixel_1 = pt2d_1[0];
                x_pixel_2 = pt2d_2[0];
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
                    iterOnObjIDMap (
                        id, 
                        y_pixel, col,
                        sight2D_list,
                        sight3D_list,
                        trID_match, 
                        trID_match_list,
                        error_list,
                        tr2d_list
                    );
                }
            }
        }
        else 
        {
            sight2D_axis.unit_vector = Pt2D(0,1);
            int y_pixel_1, y_pixel_2, min, max;

            int x_start = 0;
            int x_end = n_col;
            for (int x_pixel = x_start; x_pixel < x_end; x_pixel ++)
            {
                sight2D_axis.pt[0] = x_pixel;
                sight2D_axis.pt[1] = 0;

                // Get y_pixel (row id) from two lines 
                pt2d_1 = myMATH::crossPoint(sight2D_axis, sight2D_plus);
                pt2d_2 = myMATH::crossPoint(sight2D_axis, sight2D_minus);
                y_pixel_1 = pt2d_1[1];
                y_pixel_2 = pt2d_2[1];
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
                    iterOnObjIDMap (
                        id, 
                        row, x_pixel,
                        sight2D_list,
                        sight3D_list,
                        trID_match, 
                        trID_match_list,
                        error_list,
                        tr2d_list
                    );
                }
            }
        }

    }
    //                                  //
    // cam_cur_id = 2 ~ _check_cam_id-1 //
    //                                  //
    else
    {
        // find search region
        PixelRange search_region_all = findSearchRegion(id, sight2D_list);

        int row_start = std::min(std::max(search_region_all.row_min, 0), n_row-1);
        int row_end   = std::min(std::max(search_region_all.row_max, 1), n_row);
        int col_start = std::min(std::max(search_region_all.col_min, 0), n_col-1);
        int col_end   = std::min(std::max(search_region_all.col_max, 1), n_col);
        // std::cout << "row:" << row_start << "," << row_end << ";"
        //           << "col:" << col_start << "," << col_end
        //           << "\n";

        // iterate every pixel in the search region
        for (int i = row_start; i < row_end; i ++)
        {
            for (int j = col_start; j < col_end; j ++)
            {
                // judge whether the distances between the candidate
                // and all the lines are all within the range 
                iterOnObjIDMap (
                    id, 
                    i, j,
                    sight2D_list,
                    sight3D_list,
                    trID_match, 
                    trID_match_list,
                    error_list,
                    tr2d_list
                );
            }
        }
    }
}


void StereoMatch::iterOnObjIDMap (
    int id, 
    int row_id, int col_id,
    std::vector<Line2D> const& sight2D_list,
    std::vector<Line3D>& sight3D_list,
    std::vector<int>& trID_match, 
    std::deque<std::vector<int>>& trID_match_list,
    std::deque<double>& error_list,
    std::vector<std::vector<Tracer2D>> const& tr2d_list
)
{
    int camID_curr = _cam_list.useid_list[id];
    Pt3D pt3d;

    for (
        int k = 0; 
        k < _objID_map_list[id](row_id, col_id).size(); 
        k ++
    )
    {
        int tr_id = _objID_map_list[id](row_id, col_id)[k];

        if (tr_id == -1)
        {
            break;
        }

        bool in_range = true;
        for (int m = 0; m < id; m ++)
        {
            double dist = myMATH::distance(tr2d_list[id][tr_id]._pt_center, sight2D_list[m]);
            if (dist > _param.tor_2d)
            {
                in_range = false;
                break;
            }
        } 

        // reproject onto the previous cam, then is within error line
        if (in_range && checkReProject(id, tr_id, trID_match, tr2d_list))
        {
            trID_match.push_back(tr_id);

            // if there is still other cameras to check
            if (id < _param.check_id - 1) 
            {
                int next_id = id + 1;

                // move onto the next camera and search candidates
                findTracerMatch (
                    next_id,
                    trID_match, 
                    trID_match_list,
                    error_list,
                    tr2d_list                  
                );

                trID_match.pop_back(); // clear the match for the next candidate.
            }
            // if the current camera is the last one, then finalize the match.
            else
            {
                // test 3d distance by triangulation  
                sight3D_list[id] = _cam_list.cam_list[camID_curr].lineOfSight(tr2d_list[id][tr_id]._pt_center);
                double error_3d = 0.0;
                myMATH::triangulation(pt3d, error_3d, sight3D_list);

                if (error_3d > _param.tor_3d)
                {
                    trID_match.pop_back();
                }
                else 
                {
                    if (_param.check_id < _n_cam_use)
                    {
                        int next_id = id + 1;

                        checkTracerMatch (
                            next_id,
                            trID_match_list,
                            trID_match, 
                            pt3d,
                            error_list
                        );

                        trID_match.pop_back();
                    }
                    else 
                    {
                        trID_match_list.push_back(trID_match);

                        error_list.push_back(error_3d);

                        trID_match.pop_back();
                    }
                }
            }
        }

    }  
}


PixelRange StereoMatch::findSearchRegion (int id, std::vector<Line2D> const& sight2D_list)
{
    PixelRange search_region;
    int search_region_state = 0;

    if (sight2D_list.size() == 1)
    {
        search_region.row_max = _cam_list.cam_list[_cam_list.useid_list[id]].getNRow();
        search_region.col_max = _cam_list.cam_list[_cam_list.useid_list[id]].getNCol();
        return search_region;
    }

    // Find all crossing points
    Pt2D pt_cross;
    Pt2D perp_unit_curr;
    Pt2D perp_unit_test;
    std::vector<Line2D> line_curr_list(2);
    std::vector<Line2D> line_test_list(2);

    for (int i = 0; i < id-1; i ++)
    {
        perp_unit_curr[0] = sight2D_list[i].unit_vector[1];
        perp_unit_curr[1] = - sight2D_list[i].unit_vector[0];

        for (int j = i+1; j < id; j ++)
        {
            perp_unit_test[0] = sight2D_list[j].unit_vector[1];
            perp_unit_test[1] = - sight2D_list[j].unit_vector[0];

            // find the bounded lines for each crossing line
            line_curr_list[0].pt = sight2D_list[i].pt + perp_unit_curr * _param.tor_2d;
            line_curr_list[0].unit_vector = sight2D_list[i].unit_vector;
            line_curr_list[1].pt = sight2D_list[i].pt - perp_unit_curr * _param.tor_2d;
            line_curr_list[1].unit_vector = sight2D_list[i].unit_vector;

            line_test_list[0].pt = sight2D_list[j].pt + perp_unit_test * _param.tor_2d;
            line_test_list[0].unit_vector = sight2D_list[j].unit_vector;
            line_test_list[1].pt = sight2D_list[j].pt - perp_unit_test * _param.tor_2d;
            line_test_list[1].unit_vector = sight2D_list[j].unit_vector;

            // find the crossing points
            for (int k = 0; k < 2; k ++)
            {
                for (int l = 0; l < 2; l ++)
                {
                    pt_cross = myMATH::crossPoint(line_curr_list[k], line_test_list[l]);

                    if (search_region_state == 0)
                    {
                        search_region.row_min = int(pt_cross[1]);
                        search_region.row_max = search_region.row_min + 1;
                        search_region.col_min = int(pt_cross[0]);
                        search_region.col_max = search_region.col_min + 1;

                        search_region_state = 1;
                    }
                    else 
                    {
                        search_region.setRange(pt_cross[1], pt_cross[0]);
                    }
                }
            }
        }
    }

    return search_region;
}


// PixelRange StereoMatch::findSearchRegionOrig (int id, std::vector<Line2D> const& sight2D_list)
// {
//     Camera cam = _cam_list[cam_id];
//     PixelRange search_region;

//     if (pt_list.size() == 1)
//     {
//         // cam_id = 1, 2nd cam
//         search_region._row_max = cam.GetNpixh();
//         search_region._col_max = cam.GetNpixw();
//         return search_region;
//     }
//     else if (pt_list.size() == 2) 
//     {
//         // cam_id = 2, 3nd cam
//         Matrix<double> pt_cross (
//             FindCrossPoint (
//                 pt_list[0], 
//                 unit_list[0], 
//                 pt_list[1], 
//                 unit_list[1]
//             )
//         );

//         double cos = std::max(std::fabs(unit_list[0].Dot(unit_list[1])), MAGSMALLNUMBER);
//         double sin = std::max(std::sqrt(1-cos*cos), MAGSMALLNUMBER);

//         double width = 1.0;
//         if (cos < sin)
//         {
//             width = _tor_2d/cos;
//         }
//         else 
//         {
//             width = _tor_2d/sin;
//         }
//         width = std::min(width, 1.2*_fov_length_mm);

//         Matrix<double> pt_upright(pt_cross);
//         pt_upright(0,0) = pt_upright(0,0) + width;
//         pt_upright(1,0) = pt_upright(1,0) + width;
//         pt_upright = cam.ImgMMToPixel(pt_upright);

//         Matrix<double> pt_lowleft(pt_cross);
//         pt_lowleft(0,0) = pt_lowleft(0,0) - width;
//         pt_lowleft(1,0) = pt_lowleft(1,0) - width;
//         pt_lowleft = cam.ImgMMToPixel(pt_lowleft);

//         search_region._row_min = std::floor(pt_upright(1,0));
//         search_region._row_max = std::ceil(pt_lowleft(1,0))+1;
//         search_region._col_min = std::floor(pt_lowleft(0,0));
//         search_region._col_max = std::ceil(pt_upright(0,0))+1;

//     }
//     else if (pt_list.size() > 2)
//     {
//         // cam_id = 3, 4nd cam
//         int n = pt_list.size();

//         Matrix<double> pt_cross(3,1);
//         Matrix<double> pt_upright(3,1);
//         Matrix<double> pt_lowleft(3,1);
//         double cos = 0;
//         double sin = 0;
//         int row_min = 0;
//         int row_max = 0;
//         int col_min = 0;
//         int col_max = 0;
//         bool is_initial = false;
//         double width = 1.0;

//         for (int i = 0; i < n-1; i ++)
//         {
//             for (int j = i+1; j < n; j ++)
//             {
//                 try
//                 {
//                     pt_cross = FindCrossPoint (
//                         pt_list[i], 
//                         unit_list[i], 
//                         pt_list[j], 
//                         unit_list[j]
//                     );
//                 }
//                 catch(ErrorTypeID error)
//                 {
//                     if (error == error_parallel)
//                     {
//                         std::cout << cam_id << std::endl;
//                         continue;
//                     }
//                 }

//                 cos = std::max(std::fabs(unit_list[i].Dot(unit_list[j])), MAGSMALLNUMBER);
//                 sin = std::max(std::sqrt(1-cos*cos), MAGSMALLNUMBER);
//                 if (cos < sin)
//                 {
//                     width = _tor_2d/cos;
//                 }
//                 else 
//                 {
//                     width = _tor_2d/sin;
//                 }
//                 // width = std::min(width, 1.2*_fov_length_mm);

//                 pt_upright(0,0) = pt_upright(0,0) + width;
//                 pt_upright(1,0) = pt_upright(1,0) + width;
//                 pt_upright = cam.ImgMMToPixel(pt_upright);

//                 pt_lowleft(0,0) = pt_lowleft(0,0) - width;
//                 pt_lowleft(1,0) = pt_lowleft(1,0) - width;
//                 pt_lowleft = cam.ImgMMToPixel(pt_lowleft);

//                 row_min = std::floor(pt_upright(1,0));
//                 row_max = std::ceil(pt_lowleft(1,0))+1;
//                 col_min = std::floor(pt_lowleft(0,0));
//                 col_max = std::ceil(pt_upright(0,0))+1;

//                 if (is_initial)
//                 {
//                     search_region.SetRange(row_min, col_min);
//                     search_region.SetRange(row_max, col_max);
//                 }
//                 else
//                 {
//                     is_initial = true;
//                     search_region._row_min = row_min;
//                     search_region._row_max = row_max;
//                     search_region._col_min = col_min;
//                     search_region._col_max = col_max;
//                 }
//             }
//         }
    
//         if (! is_initial)
//         {
//             search_region._row_min = 0;
//             search_region._row_max = cam.GetNpixh();
//             search_region._col_min = 0;
//             search_region._col_max = cam.GetNpixw();
//         }
//     }

//     return search_region;
// }


bool StereoMatch::checkReProject (
    int id, 
    int tr_id,
    std::vector<int> const& tracer_id_match, 
    std::vector<std::vector<Tracer2D>> const& tr2d_list
)
{
    Line3D sight3D = _cam_list.cam_list[id].lineOfSight(tr2d_list[id][tr_id]._pt_center);

    Line2D sight2D;
    Pt2D pt2d_1;
    Pt2D pt2d_2;
    double dist = 0;

    for (int i = 0; i < id; i ++)
    {
        int cam_id = _cam_list.useid_list[i];

        pt2d_1 = _cam_list.cam_list[cam_id].project(sight3D.pt);
        pt2d_2 = _cam_list.cam_list[cam_id].project(sight3D.pt + sight3D.unit_vector);
        sight2D.pt = pt2d_1;
        sight2D.unit_vector = myMATH::createUnitVector(pt2d_1, pt2d_2);

        dist = myMATH::distance(tr2d_list[i][tracer_id_match[i]]._pt_center, sight2D);

        if (dist > _param.tor_2d)
        {
            return false;
        }
    }

    return true;
}


void StereoMatch::checkTracerMatch(
        int id, 
        std::deque<std::vector<int>>& trID_match_list,
        std::vector<int>& trID_match, 
        Pt3D& pt3d,
        std::deque<double>& error_list
    )
{

}

#endif