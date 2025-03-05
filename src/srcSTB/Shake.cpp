#ifndef SHAKE_CPP
#define SHAKE_CPP

#include "Shake.h"

void Shake::runShake(std::vector<Tracer3D>& tr3d_list, OTF const& otf, std::vector<Image> const& imgOrig_list, bool tri_only)
{
    shakeTracers(tr3d_list, otf, imgOrig_list, tri_only);
}


void Shake::checkReaptedObj(std::vector<Tracer3D> const& tr3d_list, double tol_3d)
{
    int n_tr3d = tr3d_list.size();
    _is_repeated.resize(n_tr3d, 0);
    _n_repeated = 0;

    // Remove repeated tracks
    if (_n_thread > 0)
    {
        omp_set_num_threads(_n_thread);
    }

    double repeat_thres_2 = tol_3d * tol_3d;
    for (int i = 0; i < n_tr3d-1; i ++)
    {
        if (_is_repeated[i])
        {
            continue;
        }

        #pragma omp parallel for
        for (int j = i+1; j < n_tr3d; j ++)
        {
            if (myMATH::dist2(tr3d_list[i]._pt_center, tr3d_list[j]._pt_center) < repeat_thres_2)
            {
                _is_repeated[j] = 1;
                _n_repeated ++;
            }
        }
    }
}


void Shake::shakeTracers(std::vector<Tracer3D>& tr3d_list, OTF const& otf, std::vector<Image> const& imgOrig_list, bool tri_only)
{
    // update tr2d position
    int n_tr3d = tr3d_list.size();

    if (_n_thread > 0)
    {
        omp_set_num_threads(_n_thread);
    }
    #pragma omp parallel for
    for (int i = 0; i < n_tr3d; i ++)
    {
        tr3d_list[i].projectObject2D(_cam_list.useid_list, _cam_list.cam_list);
    }

    // Initialize lists
    _imgRes_list.clear();
    _is_ghost.resize(n_tr3d, 0);
    _n_ghost = 0;

    // Initialize residue image
    int cam_id;
    for (int id = 0; id < _n_cam_use; id ++)
    {
        cam_id = _cam_list.useid_list[id];
        _imgRes_list.push_back(imgOrig_list[cam_id]);
    }

    // if only do triangulation, then skip the following steps
    if (tri_only)
    {
        calResImg(tr3d_list, otf, imgOrig_list);
        absResImg();
        return;
    }

    // Initialize score list
    _score_list.resize(n_tr3d);
    std::fill(_score_list.begin(), _score_list.end(), 1);

    double delta;
    std::vector<int> is_ignore(n_tr3d, 0);
    for (int loop = 0; loop < _n_loop; loop ++)
    {
        // update residue img
        calResImg(tr3d_list, otf, imgOrig_list);

        // update shake width
        if (loop < 1)
        {
            delta = _shake_width;
        }
        // else if (loop < 32)
        else if (loop < 5)
        {
            delta = _shake_width / std::pow(2, loop-1);  
            // delta = _shake_width / std::pow(10, loop);   
        }
        else 
        {
            delta = _shake_width / 20;
        }

        // shake each tracer
        if (_n_thread != 0)
        {
            omp_set_num_threads(_n_thread);
        }
        #pragma omp parallel
        {
            #pragma omp for
            for (int i = 0; i < n_tr3d; i ++)
            {
                if (!is_ignore[i])
                {
                    _score_list[i] = shakeOneTracer(tr3d_list[i], otf, delta, _score_list[i]);
                    // _score_list[i] = shakeOneTracerGrad(tr3d_list[i], otf, delta, _score_list[i], delta*0.1);
                }
            }
        }
    }

    // remove ghost tracers
    findGhost(tr3d_list);
    // removeGhostResidue(tr3d_list);

    // update residue image
    calResImg(tr3d_list, otf, imgOrig_list);

    // calculate absolute intensity value
    absResImg();
}


void Shake::calResImg(std::vector<Tracer3D> const& tr3d_list, OTF const& otf, std::vector<Image> const& imgOrig_list)
{
    int n_tr3d = tr3d_list.size();
    int cam_id, row_min, row_max, col_min, col_max;
    double residue;

    // initialize residue image
    for (int id = 0; id < _n_cam_use; id ++)
    {
        cam_id = _cam_list.useid_list[id];
        _imgRes_list[id] = imgOrig_list[cam_id];
    }

    // remove all tracers
    double ratio_region = 1;
    for (int i = 0; i < n_tr3d; i ++)
    {
        if (_is_ghost[i])
        {
            continue;
        }

        for (int id = 0; id < _n_cam_use; id ++)
        {
            cam_id = _cam_list.useid_list[id];

            std::vector<double> otf_param = otf.getOTFParam(cam_id, tr3d_list[i]._pt_center);

            PixelRange res_region = findRegion(
                id, 
                tr3d_list[i]._tr2d_list[id]._pt_center[1], // row_id
                tr3d_list[i]._tr2d_list[id]._pt_center[0], // col_id
                tr3d_list[i]._tr2d_list[id]._r_px * ratio_region
            );
        
            row_min = res_region.row_min;
            row_max = res_region.row_max;
            col_min = res_region.col_min;
            col_max = res_region.col_max;

            for (int row = row_min; row < row_max; row ++)
            {
                for (int col = col_min; col < col_max; col ++)
                {
                    residue = imgOrig_list[cam_id](row, col) - gaussIntensity(col, row, tr3d_list[i]._tr2d_list[id]._pt_center, otf_param);
                    
                    if (!std::isfinite(residue))
                    {
                        std::cerr << "Shake::calResImg error at line " << __LINE__ << ".\n";
                        std::cerr << "id = " << id << ", cam_id = " << cam_id << "\n";
                        std::cerr << "(row,col) = " << "(" << row << "," << col << ")\n";
                        std::cerr << "residue = " << residue << "\n";

                        std::cout << "otf_param = ";
                        for (int i = 0; i < otf_param.size(); i++)
                        {
                            std::cerr << otf_param[i] << ',';
                        }
                        std::cerr << std::endl;
                        
                        std::cerr << "pt2d = " << std::endl;
                        tr3d_list[i]._tr2d_list[id]._pt_center.print();

                        std::cerr << "pt3d = " << std::endl;
                        tr3d_list[i]._pt_center.print();

                        throw error_range;
                    }

                    // choose the min residue as the value
                    // TODO: try other ways
                    // if (residue < _imgRes_list[id](row, col))
                    // {
                    //     _imgRes_list[id](row, col) = residue;
                    // }
                    _imgRes_list[id](row, col) = residue;
                }
            }
        }
    }
}


PixelRange Shake::findRegion(int id, double row, double col, double half_width_px)
{
    int cam_id = _cam_list.useid_list[id];
    int n_row = _cam_list.cam_list[cam_id].getNRow();
    int n_col = _cam_list.cam_list[cam_id].getNCol();

    PixelRange region;

    // region.row_min = std::min(std::max(0, row-half_width_px), n_row-1);

    // region.col_min = std::min(std::max(0, col-half_width_px), n_col-1);

    // region.row_max = std::max(std::min(n_row, row+half_width_px+1), 1);

    // region.col_max = std::max(std::min(n_col, col+half_width_px+1), 1);

    region.row_min = std::max(0, int(std::floor(row-half_width_px)));

    region.col_min = std::max(0, int(std::floor(col-half_width_px)));

    region.row_max = std::min(n_row, int(std::ceil(row+half_width_px+1)));

    region.col_max = std::min(n_col, int(std::ceil(col+half_width_px+1)));

    return region;
}


double Shake::gaussIntensity(int x, int y, Pt2D const& pt2d, std::vector<double> const& otf_param)
{
    double dx = x - pt2d[0];
    double dy = y - pt2d[1];
    double xx =  dx * std::cos(otf_param[3]) + dy * std::sin(otf_param[3]);
    double yy = -dx * std::sin(otf_param[3]) + dy * std::cos(otf_param[3]);
    double value = otf_param[0] * std::exp(- otf_param[1] * (xx*xx) - otf_param[2] * (yy*yy));
    return std::max(0.0, value);
}


double Shake::shakeOneTracer(Tracer3D& tr3d, OTF const& otf, double delta, double score_old)
{
    ImgAugList imgAug_list; // augmented image list
    
    // calculate augmented image for each camera 
    // pixel range radius: _r_px
    int ratio_region = 2;
    int cam_id;
    for (int id = 0; id < _n_cam_use; id ++)
    {
        cam_id = _cam_list.useid_list[id];

        // Get the augmented image range   
        PixelRange region = findRegion(
            id, 
            tr3d._tr2d_list[id]._pt_center[1], // row
            tr3d._tr2d_list[id]._pt_center[0], // col
            tr3d._tr2d_list[id]._r_px * ratio_region
        );
        imgAug_list.region_list.push_back(region);

        // Get the range to calculate the projection intensity
        PixelRange int_region = findRegion(
            id, 
            tr3d._tr2d_list[id]._pt_center[1], // row
            tr3d._tr2d_list[id]._pt_center[0], // col
            tr3d._tr2d_list[id]._r_px
        );
        

        // create a particle reproj image (I_p) matrix in the pixel range
        int n_row = region.getNumOfRow();
        int n_col = region.getNumOfCol();
        if (n_row <= 0 || n_col <= 0)
        {
            n_row = 1;
            n_col = 1;
        }
        Image aug_img(
            n_row,
            n_col, 
            0
        );

        std::vector<double> otf_para = otf.getOTFParam(cam_id, tr3d._pt_center);

        int i = 0;
        for (int row = region.row_min; row < region.row_max; row ++)
        {
            int j = 0;
            for (int col = region.col_min; col < region.col_max; col ++)
            {
                double value = 0;
                
                bool judge = row >= int_region.row_min && 
                             row < int_region.row_max && 
                             col >= int_region.col_min && 
                             col < int_region.col_max;
                if (judge)
                {
                    value = gaussIntensity(
                        col, row,
                        tr3d._tr2d_list[id]._pt_center,
                        otf_para
                    );
                }

                // Creating a particle augmented residual image: (res+p)
                aug_img(i, j) = std::max( 
                    0.0,
                    std::min(double(_cam_list.intensity_max[cam_id]), 
                             _imgRes_list[id](row,col) + value)
                );
                j++;
            }
            i++;
        }

        imgAug_list.img_list.push_back(aug_img);
    }

    // Update the particle position, imgAug and search range
    double residue = updateTracer(tr3d, imgAug_list, otf, delta);

    // update the tracer score
    // sum up the score of all cam 
    // to delete ghost tracer 
    double tr_score = calTracerScore(
        tr3d,
        imgAug_list, 
        otf,
        score_old
    );
    return tr_score;

    // return residue;
}


double Shake::shakeOneTracerGrad(Tracer3D& tr3d, OTF const& otf, double delta, double score_old, double lr)
{
    ImgAugList imgAug_list; // augmented image list
    
    int cam_id;
    double ratio_region = 2;
    for (int id = 0; id < _n_cam_use; id ++)
    {
        cam_id = _cam_list.useid_list[id];

        PixelRange region = findRegion(
            id, 
            tr3d._tr2d_list[id]._pt_center[1], // row
            tr3d._tr2d_list[id]._pt_center[0], // col
            tr3d._tr2d_list[id]._r_px * ratio_region
        );

        // create a particle reproj image (I_p) matrix in the pixel range
        // Image aug_img(
        //     region.getNumOfRow(),
        //     region.getNumOfCol(), 
        //     0
        // );
        Image aug_img;
        if (region.getNumOfRow() > 0 && region.getNumOfCol() > 0)
        {
            aug_img = Image(
                region.getNumOfRow(),
                region.getNumOfCol(), 
                0
            );
        }
        
        std::vector<double> otf_para = otf.getOTFParam(cam_id, tr3d._pt_center);

        int i = 0;
        for (int row = region.row_min; row < region.row_max; row ++)
        {
            int j = 0;
            for (int col = region.col_min; col < region.col_max; col ++)
            {
                double value = gaussIntensity(
                    col, row,
                    tr3d._tr2d_list[id]._pt_center,
                    otf_para
                );

                // Creating a particle augmented residual image: (res+p)
                aug_img(i, j) = std::max( 
                    0.0,
                    std::min(double(_cam_list.intensity_max[cam_id]), 
                             _imgRes_list[id](row,col) + value)
                );
                j++;
            }
            i++;
        }

        imgAug_list.region_list.push_back(region);
        imgAug_list.img_list.push_back(aug_img);
    }

    // Update the particle position, imgAug and search range
    double residue = updateTracerGrad(tr3d, imgAug_list, otf, delta, lr);

    // update the tracer score
    // sum up the score of all cam 
    // to delete ghost tracer 
    double tr_score = calTracerScore(
        tr3d,
        imgAug_list, 
        otf,
        score_old
    );
    return tr_score;

    // return residue;
}


double Shake::updateTracer(Tracer3D& tr3d, ImgAugList& imgAug_list, OTF const& otf, double delta)
{
    std::vector<double> delta_list(3);
    delta_list[0] = - delta;
    delta_list[1] = 0;
    delta_list[2] = delta;

    std::vector<double> array_list(4, 0);
    std::vector<double> array_list_fit(3, 0);
    std::vector<double> coeff(3, 0);
    std::vector<double> residue_list(4, 0);
    std::vector<double> residue_list_fit(3, 0);

    // shaking on x,y,z direction
    int min_id;
    double residue = 0.0;
    Tracer3D tr3d_temp(tr3d);

    for (int i = 0; i < 3; i ++)
    {
        for (int j = 0; j < 3; j ++)
        {
            array_list[j] = tr3d._pt_center[i] + delta_list[j];
            array_list_fit[j] = array_list[j];

            tr3d_temp._pt_center[i] = array_list[j];
            
            // update tr3d 2d match
            tr3d_temp.projectObject2D(_cam_list.useid_list, _cam_list.cam_list);

            // update imgAug_list and region_list
            // updateImgAugList(imgAug_list, tr3d_temp);   

            residue_list[j] = calPointResidue(
                tr3d_temp,
                imgAug_list,
                otf
            );
            residue_list_fit[j] = residue_list[j];
        }
        
        // residue = coeff[0] + coeff[1] * x + coeff[2] * x^2
        myMATH::polyfit(coeff, array_list_fit, residue_list_fit, 2);
        // coeff[2] = (((residue_list[1] - residue_list[0]) / (array_list[1] - array_list[0])) - ((residue_list[2] - residue_list[0]) / (array_list[2] - array_list[0]))) / (array_list[1] - array_list[2]); //a
        // coeff[1] = ((residue_list[1] - residue_list[0]) / (array_list[1] - array_list[0])) - coeff[2]*(array_list[1] + array_list[0]); //b
        // coeff[0] = residue_list[2] - coeff[2] * std::pow(array_list[2], 2) - coeff[1]*array_list[2]; //c


        if (coeff[2] == 0)
        {
            if (coeff[1] == 0)
            {
                min_id = 1;
            }
            else if (coeff[1] > 0)
            {
                min_id = 0;
            }
            else
            {
                min_id = 2;
            }
        }
        else 
        {
            array_list[3] = - coeff[1] / (2 * coeff[2]);
            if (array_list[3]>array_list[0] && array_list[3]<array_list[2])
            {
                tr3d_temp._pt_center[i] = array_list[3];

                // update tr3d 2d match
                tr3d_temp.projectObject2D(_cam_list.useid_list, _cam_list.cam_list);
                
                // update imgAug_list and region_list
                // updateImgAugList(imgAug_list, tr3d_temp);

                residue_list[3] = calPointResidue(
                    tr3d_temp,
                    imgAug_list,
                    otf
                );
            }
            else
            {
                residue_list[3] = residue_list[0] + residue_list[1] + residue_list[2] + 1; // set a maximum residue value
            }

            min_id = std::min_element(residue_list.begin(), residue_list.end()) - residue_list.begin();
        }
        
        // update tr3d 2d match
        tr3d._pt_center[i] = array_list[min_id];
        residue = residue_list[min_id];
        // tr3d.projectObject2D(_cam_list.useid_list, _cam_list.cam_list);
        // // update imgAug_list and region_list
        // updateImgAugList(imgAug_list, tr3d);
    }

    tr3d.projectObject2D(_cam_list.useid_list, _cam_list.cam_list);
    // updateImgAugList(imgAug_list, tr3d);

    return residue;
}


double Shake::updateTracerGrad(Tracer3D& tr3d, ImgAugList& imgAug_list, OTF const& otf, double delta, double lr)
{
    // AxisLimit limit(
    //     tr3d._pt_center[0]-delta, tr3d._pt_center[0]+delta,
    //     tr3d._pt_center[1]-delta, tr3d._pt_center[1]+delta,
    //     tr3d._pt_center[2]-delta, tr3d._pt_center[2]+delta
    // );

    double residue_old, residue_new;
    Pt3D shift3d, grad;

    // initialize shift
    std::mt19937 gen(1);
    std::bernoulli_distribution distribution(0.5);
    shift3d[0] = int(distribution(gen)) * 2 - 1;
    shift3d[1] = int(distribution(gen)) * 2 - 1;
    shift3d[2] = int(distribution(gen)) * 2 - 1;
    shift3d = shift3d / shift3d.norm() * delta;

    // calculate current residue
    tr3d.projectObject2D(_cam_list.useid_list, _cam_list.cam_list);
    residue_old = calPointResidue(tr3d, imgAug_list, otf);

    double residue_orig = residue_old; // for debug
    Pt3D pt_center_orig = tr3d._pt_center; // for debug

    double tor = 1e-6;
    int record_step = 0;
    double norm_grad = std::pow(40/40, 2);
    for (int step = 0; step < 50000; step ++)
    {
        // std::cout << std::endl;
        // std::cout << "step = " << step << ", residue_old = " << residue_old << std::endl;
        // std::cout << "shift3d = " << shift3d[0] << ", " << shift3d[1] << ", " << shift3d[2] << std::endl;
        // std::cout << "tr3d._pt_center = " << tr3d._pt_center[0] << ", " << tr3d._pt_center[1] << ", " << tr3d._pt_center[2] << std::endl;

        // update tr3d
        // int check_step = 0;
        // while (
        //     ! limit.check(
        //         tr3d._pt_center[0] + shift3d[0],
        //         tr3d._pt_center[1] + shift3d[1],
        //         tr3d._pt_center[2] + shift3d[2]
        //     )
        //     && check_step < 10
        // )
        // {
        //     shift3d /= 10;
        //     check_step ++;
        // }
        tr3d._pt_center += shift3d;

        // update tr2d
        tr3d.projectObject2D(_cam_list.useid_list, _cam_list.cam_list);

        // update imgAug_list and region_list
        // updateImgAugList(imgAug_list, tr3d);

        // calculate new residue
        residue_new = calPointResidue(tr3d, imgAug_list, otf);

        // calculate grad
        for (int i = 0; i < 3; i ++)
        {
            grad[i] = (residue_new - residue_old) / shift3d[i];
        }

        // update shift 3d
        norm_grad = 1 / grad.norm();
        shift3d = grad * norm_grad * (-lr);

        // check if the possition is converged
        if (shift3d.norm() < tor)
        {
            break;
        }

        // update residue_old
        residue_old = residue_new;
        record_step ++;
    }

    
    if (residue_new < residue_orig)
    {
        // std::cout << "step = " << record_step << ", shift3d.mag = " << shift3d.norm() << ", residue_orig = " << residue_orig << ", residue_new = " << residue_new <<  std::endl;
    }
    else
    {
        tr3d._pt_center = pt_center_orig;
        residue_new = residue_orig;
    }
    std::cout << "step = " << record_step << ", shift3d.mag = " << shift3d.norm() << ", residue_orig = " << residue_orig << ", residue_new = " << residue_new <<  std::endl;
    
    return residue_new;
}


void Shake::updateImgAugList (ImgAugList& imgAug_list, Tracer3D const& tr3d)
{
    double ratio_region = 2;

    // update Augemented image and search range
    for (int id = 0; id < _n_cam_use; id ++)
    {
        // Get the new range based on the updated 3d position
        PixelRange region = findRegion(
            id, 
            tr3d._tr2d_list[id]._pt_center[1], 
            tr3d._tr2d_list[id]._pt_center[0], 
            tr3d._tr2d_list[id]._r_px * ratio_region
        );
        int n_col_new = region.getNumOfCol();
        int n_row_new = region.getNumOfRow();
        int n_col_old = imgAug_list.region_list[id].getNumOfCol();
        int n_row_old = imgAug_list.region_list[id].getNumOfRow();

        // save the orig imgAug and update imgAug based on new range
        Image imgAug_old = imgAug_list.img_list[id];
        if (n_col_new != n_col_old || 
            n_row_new != n_row_old)
        {
            if (n_col_new <= 0 || n_row_new <= 0)
            {
                n_col_new = 1;
                n_row_new = 1;
            }

            imgAug_list.img_list[id] = Image(n_row_new, n_col_new, 0);
        }
        
        // record the old range
        int row_min_old = imgAug_list.region_list[id].row_min;
        int row_max_old = imgAug_list.region_list[id].row_max;
        int col_min_old = imgAug_list.region_list[id].col_min;
        int col_max_old = imgAug_list.region_list[id].col_max;

        // update search range
        imgAug_list.region_list[id] = region;

        // update imgAug based on new range
        int i = 0; 
        for (int row = region.row_min; row < region.row_max; row ++)
        {
            int j = 0; 
            for (int col = region.col_min; col < region.col_max; col ++)
            {
                bool judge = row <  row_max_old && 
                             col <  col_max_old &&
                             row >= row_min_old &&
                             col >= col_min_old;
                if (judge)
                {
                    imgAug_list.img_list[id](i,j) = imgAug_old(row-row_min_old, col-col_min_old);
                }
                else
                {
                    imgAug_list.img_list[id](i,j) = std::max(0.0, _imgRes_list[id](row, col));
                }
                j ++;
            }
            i ++;
        }   
    }    
}


// Try with image cross correlation
double Shake::calPointResidue (Tracer3D const& tr3d, ImgAugList const& imgAug_list, OTF const& otf)
{
    double residue = 0;
    int cam_id;

    for (int id = 0; id < _n_cam_use; id ++)
    {
        cam_id = _cam_list.useid_list[id];

        // get otf param
        std::vector<double> otf_param = otf.getOTFParam(cam_id, tr3d._pt_center);

        // get the region for calculating projection intensity
        PixelRange int_region = findRegion(
            id, 
            tr3d._tr2d_list[id]._pt_center[1], // row
            tr3d._tr2d_list[id]._pt_center[0], // col
            tr3d._tr2d_list[id]._r_px
        );

        // Calculating the residual for updated 3D position
        //  residue = sum(sum((I_p - I_o)^2), all cam) / n_cam
        int i = 0;
        for (int row = imgAug_list.region_list[id].row_min; 
                 row < imgAug_list.region_list[id].row_max; 
                 row ++)
        {
            int j = 0;
            for (int col = imgAug_list.region_list[id].col_min; 
                     col < imgAug_list.region_list[id].col_max; 
                     col ++)
            {
                double value = 0;

                // calculate gaussian intensity when distance is r_px
                bool judge = row >= int_region.row_min && 
                             row < int_region.row_max && 
                             col >= int_region.col_min && 
                             col < int_region.col_max;
                if (judge)
                {
                    value = gaussIntensity(
                        col, row,
                        tr3d._tr2d_list[id]._pt_center,
                        otf_param
                    );
                }
                
                // value = gaussIntensity(
                //     col, row,
                //     tr3d._tr2d_list[id]._pt_center,
                //     otf_param
                // );

                residue += std::pow(
                    imgAug_list.img_list[id](i,j) - value,
                    2
                );

                j ++;
                // npts ++;
            }
            i ++;
        }
    }

    // residue /= npts;
    return residue;
}


double Shake::calTracerScore (Tracer3D const& tr3d, ImgAugList const& imgAug_list, OTF const& otf, double score)
{
    int cam_id;


    // get range for calculating score (radius = r_px)
    std::vector<PixelRange> score_region(_n_cam_use);
    for (int id = 0; id < _n_cam_use; id ++)
    {
        cam_id = _cam_list.useid_list[id];

        // Get the actual pixel range   
        score_region[id] = findRegion(
            id, 
            tr3d._tr2d_list[id]._pt_center[1], // row
            tr3d._tr2d_list[id]._pt_center[0], // col
            tr3d._tr2d_list[id]._r_px
        );
    }


    // update particle intensity 
    // (ignoring the camera with highest local peak intensity)
    // Calculate intensity for each camera: 
    //  if the tracer is out of the image, directly set the intensity as 0
    //  then calculate the new score without the outlier
    int n_outRange = 0;
    std::vector<int> is_select(_n_cam_use, 1);

    std::vector<double> peakInt_list(_n_cam_use, 0);
    for (int id = 0; id < _n_cam_use; id ++)
    {
        // if the tracer is out of the image, do not select this camera
        bool judge = score_region[id].getNumOfRow() <= 0 || 
                     score_region[id].getNumOfCol() <= 0;
        if (judge)
        {
            n_outRange ++;
            is_select[id] = 0;

            if (_n_cam_use-n_outRange < 2)
            {
                return score * 0.1;
            } 

            continue;
        }
        
        // find the peak intensity in score region
        for (int row = score_region[id].row_min; 
                 row < score_region[id].row_max; 
                 row ++)
        {
            for (int col = score_region[id].col_min; 
                     col < score_region[id].col_max; 
                     col ++)
            {
                int i = row - imgAug_list.region_list[id].row_min;
                int j = col - imgAug_list.region_list[id].col_min;

                bool judge = i >= 0 && 
                             i < imgAug_list.region_list[id].getNumOfRow() &&
                             j >= 0 && 
                             j < imgAug_list.region_list[id].getNumOfCol();
                
                if (judge)
                {
                    peakInt_list[id] = std::max(peakInt_list[id], imgAug_list.img_list[id](i, j));
                }
                else
                {
                    peakInt_list[id] = std::max(peakInt_list[id], _imgRes_list[id](row, col));
                }
            }
        }
    }

    // ignore the camera with highest local peak intensity
    if (_n_cam_use - n_outRange > 2)
    {
        int cam_id_ignore = std::max_element(peakInt_list.begin(), peakInt_list.end()) - peakInt_list.begin();
        is_select[cam_id_ignore] = 0;
    }


    // update score
    std::vector<double> numerator_list(_n_cam_use, 0);
    std::vector<double> denominator_list(_n_cam_use, 0);
    for (int id = 0; id < _n_cam_use; id ++)
    {
        if (!is_select[id])
        {
            continue;
        }

        std::vector<double> otf_param = otf.getOTFParam(
            _cam_list.useid_list[id],
            tr3d._pt_center
        );

        for (int row = score_region[id].row_min; 
                 row < score_region[id].row_max; 
                 row ++)
        {
            for (int col = score_region[id].col_min; 
                     col < score_region[id].col_max; 
                     col ++)
            {
                int i = row - imgAug_list.region_list[id].row_min;
                int j = col - imgAug_list.region_list[id].col_min;

                bool judge = i >= 0 && 
                             i < imgAug_list.region_list[id].getNumOfRow() &&
                             j >= 0 && 
                             j < imgAug_list.region_list[id].getNumOfCol();

                if (judge)
                {
                    numerator_list[id] += imgAug_list.img_list[id](i, j); 
                }
                else
                {
                    numerator_list[id] += _imgRes_list[id](row, col);
                
                }
                
                denominator_list[id] += gaussIntensity(col, row, tr3d._tr2d_list[id]._pt_center, otf_param);
            }
        }
    }

    // calculate numerator and denominator
    double numerator = 0.0;     
    double denominator = 0.0;
    for (int id = 0; id < _n_cam_use; id ++)
    {
        numerator += numerator_list[id];
        denominator += denominator_list[id];   
    }    
    
    if (denominator > SMALLNUMBER)
    {
        score *= std::sqrt(std::abs(numerator/denominator));
    }
    return score;
}


void Shake::findGhost(std::vector<Tracer3D>& tr3d_list)
{
    int n_tr3d = tr3d_list.size();

    // find particles that are close to each other
    // _is_repeated.resize(n_tr3d, 0);
    checkReaptedObj(tr3d_list, _tol_3d);

    // remove outliers
    double sum = 0;
    int n_mean = 0;
    double mean = 0;
    for (int i = 0; i < n_tr3d; i ++)
    {
        if (_score_list[i] > 0 && !_is_repeated[i])
        {
            sum += _score_list[i];
            n_mean ++;
        }
    }
    n_mean = n_tr3d;
    mean = sum / n_mean;

    for (int i = 0; i < n_tr3d; i ++)
    {
        if (_score_list[i] > 30 * mean && !_is_repeated[i])
        {
            sum -= _score_list[i];
            n_mean --;
        }
    }
    mean = sum / n_mean;

    // std::cout << "Shake score mean: " << mean << std::endl;

    // remove ghost tracers if the score is less than _min_score*mean
    _n_ghost = 0;
    for (int i = 0; i < n_tr3d; i ++)
    {
        if (_score_list[i] < _score_min * mean || _is_repeated[i])
        {
            _is_ghost[i] = 1;
            _n_ghost ++;
        }
        else
        {
            _is_ghost[i] = 0;
        }
    }

    #ifdef DEBUG
    std::cout << "\tShake::findGhost: find " << _n_ghost << " ghost tracers, including " << _n_repeated << " repeated tracers." << std::endl;
    #endif
}


void Shake::absResImg ()
{
    for (int id = 0; id < _n_cam_use; id ++)
    {
        int cam_id = _cam_list.useid_list[id];
        int n_row = _cam_list.cam_list[cam_id].getNRow();
        int n_col = _cam_list.cam_list[cam_id].getNCol();

        for (int row = 0; row < n_row; row ++)
        {
            for (int col = 0; col < n_col; col ++)
            {
                if (!std::isfinite(_imgRes_list[id](row, col)))
                {
                    std::cerr << "Shake::absResImg: "
                              << "id = " << id << ", "
                              << "cam_id = " << cam_id << ", "
                              << "(row, col) = " 
                              << "(" << row << ","
                              << col << "), "
                              << "_imgRes_list[id](row, col) = "
                              << _imgRes_list[id](row, col)
                              << std::endl;
                    throw error_range;
                }

                _imgRes_list[id](row, col) = std::max(0.0, _imgRes_list[id](row, col));
            }
        }
    }
}

#endif