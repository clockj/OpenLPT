#ifndef SHAKE_CPP
#define SHAKE_CPP

#include "Shake.h"

void Shake::runShake(std::vector<Tracer3D>& tr3d_list, OTF const& otf, std::vector<Image> const& imgOrig_list, bool tri_only)
{
    shakeTracers(tr3d_list, otf, imgOrig_list, tri_only);
}


void Shake::shakeTracers(std::vector<Tracer3D>& tr3d_list, OTF const& otf, std::vector<Image> const& imgOrig_list, bool tri_only)
{
    // update tr2d position
    int n_tr3d = tr3d_list.size();
    for (int i = 0; i < n_tr3d; i ++)
    {
        tr3d_list[i].projectObject2D(_cam_list.useid_list, _cam_list.cam_list);
    }

    // if only do triangulation, then skip the following steps
    if (tri_only)
    {
        calResImg(tr3d_list, otf, imgOrig_list);
        absResImg();
        return;
    }

    // Clear lists
    _imgRes_list.clear();
    _objID_keep.clear();
    _objID_remove.clear();
    _n_ghost = 0;

    // Initialize score list
    _score_list.resize(n_tr3d);
    std::fill(_score_list.begin(), _score_list.end(), 1);

    // Initialize residue image
    int cam_id;
    for (int id = 0; id < _n_cam_use; id ++)
    {
        cam_id = _cam_list.useid_list[id];
        _imgRes_list.push_back(imgOrig_list[cam_id]);
    }

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
        else if (loop < 32)
        {
            // delta = _shake_width / std::pow(2, loop);  
            delta = _shake_width / std::pow(10, loop);   
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
                    // _score_list[i] = shakeOneTracerGrad(tr3d_list[i], otf, delta, 1);
                }
            }
        }
    }

    // remove ghost tracers
    removeGhost(tr3d_list);
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

    // remove all tracers
    for (int i = 0; i < n_tr3d; i ++)
    {
        for (int id = 0; id < _n_cam_use; id ++)
        {
            cam_id = _cam_list.useid_list[id];

            std::vector<double> otf_param = otf.getOTFParam(cam_id, tr3d_list[i]._pt_center);

            PixelRange region = findRegion(
                id, 
                tr3d_list[i]._tr2d_list[id]._pt_center[1], // row_id
                tr3d_list[i]._tr2d_list[id]._pt_center[0], // col_id
                tr3d_list[i]._tr2d_list[id]._r_px
            );
        
            row_min = region.row_min;
            row_max = region.row_max;
            col_min = region.col_min;
            col_max = region.col_max;

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
                    if (residue < _imgRes_list[id](row, col))
                    {
                        _imgRes_list[id](row, col) = residue;
                    }
                }
            }
        }
    }
}


PixelRange Shake::findRegion(int id, int row, int col, int half_width_px)
{
    int cam_id = _cam_list.useid_list[id];
    int n_row = _cam_list.cam_list[cam_id].getNRow();
    int n_col = _cam_list.cam_list[cam_id].getNCol();

    PixelRange region;

    region.row_min = std::min(std::max(0, row-half_width_px), n_row-1);

    region.col_min = std::min(std::max(0, col-half_width_px), n_col-1);

    region.row_max = std::max(std::min(n_row, row+half_width_px+1), 1);

    region.col_max = std::max(std::min(n_col, col+half_width_px+1), 1);

    return region;
}


double Shake::gaussIntensity(int x, int y, Pt2D const& pt2d, std::vector<double> const& otf_param)
{
    double dx = x - pt2d[0];
    double dy = y - pt2d[1];
    double xx =  dx * std::cos(otf_param[3]) + dy * std::sin(otf_param[3]);
    double yy = -dx * std::sin(otf_param[3]) + dy * std::cos(otf_param[3]);
    double value = otf_param[0] * std::exp(-otf_param[1] * (xx*xx) - otf_param[2] * (yy*yy));
    return std::max(0.0, value);
}


double Shake::shakeOneTracer(Tracer3D& tr3d, OTF const& otf, double delta, double score_old)
{
    std::vector<PixelRange> region_list(_n_cam_use);
    std::vector<Image> imgAug_list; // augmented image list
    
    int cam_id;
    for (int id = 0; id < _n_cam_use; id ++)
    {
        cam_id = _cam_list.useid_list[id];

        // Get the actual pixel range   
        // TODO: try 1.5 search range   
        // tr3d._tr2d_list[id]._pt_center.print(); 
         
        PixelRange region = findRegion(
            id, 
            tr3d._tr2d_list[id]._pt_center[1], // row
            tr3d._tr2d_list[id]._pt_center[0], // col
            tr3d._tr2d_list[id]._r_px
        );

        // create a particle reproj image (I_p) matrix in the pixel range
        Image aug_img(
            region.getNumOfRow(),
            region.getNumOfCol(), 
            0
        );

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

        region_list[id] = region;
        imgAug_list.push_back(aug_img);
    }

    // Update the particle position, imgAug and search range
    double residue = updateTracer(tr3d, imgAug_list, region_list, otf, delta);

    // update the tracer score
    // sum up the score of all cam 
    // to delete ghost tracer 
    double tr_score = calTracerScore(
        tr3d,
        region_list,
        imgAug_list, 
        otf,
        score_old
    );
    return tr_score;

    // return residue;
}


double Shake::shakeOneTracerGrad(Tracer3D& tr3d, OTF const& otf, double delta, double lr)
{
    std::vector<PixelRange> region_list(_n_cam_use);
    std::vector<Image> imgAug_list; // augmented image list
    
    int cam_id;
    for (int id = 0; id < _n_cam_use; id ++)
    {
        cam_id = _cam_list.useid_list[id];

        PixelRange region = findRegion(
            id, 
            tr3d._tr2d_list[id]._pt_center[1], // row
            tr3d._tr2d_list[id]._pt_center[0], // col
            tr3d._tr2d_list[id]._r_px
        );

        // create a particle reproj image (I_p) matrix in the pixel range
        Image aug_img(
            region.getNumOfRow(),
            region.getNumOfCol(), 
            0
        );

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

        region_list[id] = region;
        imgAug_list.push_back(aug_img);
    }

    // Update the particle position, imgAug and search range
    double residue = updateTracerGrad(tr3d, imgAug_list, region_list, otf, delta, lr);

    return residue;
}


double Shake::updateTracer(Tracer3D& tr3d, std::vector<Image>& imgAug_list, std::vector<PixelRange>& region_list, OTF const& otf, double delta)
{
    std::vector<double> delta_list(3);
    delta_list[0] = - delta;
    delta_list[1] = 0;
    delta_list[2] = delta;

    std::vector<double> array_list(4, 0);
    std::vector<double> coeff(3, 0);
    std::vector<double> residue_list(4, 0);

    // shaking on x,y,z direction
    double limit_min, limit_max;
    int min_id;
    double residue = 0.0;
    for (int i = 0; i < 3; i ++)
    {
        switch (i)
        {
            case 0:
                limit_min = otf._param.boundary.x_min;
                limit_max = otf._param.boundary.x_max;
                break;
            case 1:
                limit_min = otf._param.boundary.y_min;
                limit_max = otf._param.boundary.y_max;
                break;
            case 2:
                limit_min = otf._param.boundary.z_min;
                limit_max = otf._param.boundary.z_max;
                break;
            default:
                std::cerr << "Shake::updateTracer: error at line " << __LINE__ << ".\n";
                throw error_range;
        }

        for (int j = 0; j < 3; j ++)
        {
            array_list[j] = tr3d._pt_center[i] + delta_list[j];
            if (array_list[j] < limit_min)
            {
                array_list[j] = limit_min;
            }
            else if (array_list[j] > limit_max)
            {
                array_list[j] = limit_max;
            }

            tr3d._pt_center[i] = array_list[j];
            
            // // update tr3d 2d match
            // tr3d.projectObject2D(_cam_list.useid_list, _cam_list.cam_list);
            // // update imgAug_list and region_list
            // updateImgAugList(imgAug_list, region_list, tr3d);   

            residue_list[j] = calPointResidue(
                tr3d._pt_center,
                region_list,
                imgAug_list,
                otf
            );
        }
        
        // residue = coeff[0] + coeff[1] * x + coeff[2] * x^2
        myMATH::polyfit(coeff, array_list, residue_list, 2);

        if (std::fabs(coeff[2]) < SMALLNUMBER)
        {
            if (std::fabs(coeff[1]) < SMALLNUMBER)
            {
                tr3d._pt_center[i] = array_list[1];
                residue = residue_list[1];
            }
            else if (coeff[1] > 0)
            {
                tr3d._pt_center[i] = array_list[0];
                residue = residue_list[0];
            }
            else
            {
                tr3d._pt_center[i] = array_list[2];
                residue = residue_list[2];
            }
        }
        else 
        {
            array_list[3] = - coeff[1] / (2 * coeff[2]);
            if (array_list[3]>array_list[0] && array_list[3]<array_list[2])
            {
                tr3d._pt_center[i] = array_list[3];
                residue_list[3] = calPointResidue(
                    tr3d._pt_center,
                    region_list,
                    imgAug_list,
                    otf
                );
            }
            else
            {
                residue_list[3] = residue_list[0] + residue_list[1] + residue_list[2] + 1; // set a maximum residue value
            }

            min_id = std::min_element(residue_list.begin(), residue_list.end()) - residue_list.begin();
            tr3d._pt_center[i] = array_list[min_id];
            residue = residue_list[min_id];
        }
        
        // update tr3d 2d match
        tr3d.projectObject2D(_cam_list.useid_list, _cam_list.cam_list);
        // update imgAug_list and region_list
        updateImgAugList(imgAug_list, region_list, tr3d);
    }

    // tr3d.projectObject2D(_cam_list.useid_list, _cam_list.cam_list);
    // updateImgAugList(imgAug_list, region_list, tr3d);

    return residue;
}


double Shake::updateTracerGrad(Tracer3D& tr3d, std::vector<Image>& imgAug_list, std::vector<PixelRange>& region_list, OTF const& otf, double delta, double lr)
{
    AxisLimit limit(
        tr3d._pt_center[0]-delta, tr3d._pt_center[0]+delta,
        tr3d._pt_center[1]-delta, tr3d._pt_center[1]+delta,
        tr3d._pt_center[2]-delta, tr3d._pt_center[2]+delta
    );

    double residue_old, residue_new;
    Pt3D shift3d, grad;

    // initialize shift
    std::mt19937 gen(1);
    std::bernoulli_distribution distribution(0.5);
    shift3d[0] = int(distribution(gen)) * 2 - 1;
    shift3d[1] = int(distribution(gen)) * 2 - 1;
    shift3d[2] = int(distribution(gen)) * 2 - 1;
    shift3d = shift3d / shift3d.norm() * delta * 0.1;

    // calculate current residue
    residue_old = calPointResidue(tr3d._pt_center, region_list, imgAug_list, otf);

    double residue_orig = residue_old; // for debug
    Pt3D pt_center_orig = tr3d._pt_center; // for debug

    double tor = 1e-4;
    int record_step = 0;
    double norm_grad = std::pow(40/40, 2);
    for (int step = 0; step < 100; step ++)
    {
        // std::cout << std::endl;
        // std::cout << "step = " << step << ", residue_old = " << residue_old << std::endl;
        // std::cout << "shift3d = " << shift3d[0] << ", " << shift3d[1] << ", " << shift3d[2] << std::endl;
        // std::cout << "tr3d._pt_center = " << tr3d._pt_center[0] << ", " << tr3d._pt_center[1] << ", " << tr3d._pt_center[2] << std::endl;

        // update tr3d
        int check_step = 0;
        while (
            ! limit.check(
                tr3d._pt_center[0] + shift3d[0],
                tr3d._pt_center[1] + shift3d[1],
                tr3d._pt_center[2] + shift3d[2]
            )
            && check_step < 10
        )
        {
            shift3d /= 10;
            check_step ++;
        }
        tr3d._pt_center += shift3d;

        // update tr2d
        tr3d.projectObject2D(_cam_list.useid_list, _cam_list.cam_list);

        // update imgAug_list and region_list
        updateImgAugList(imgAug_list, region_list, tr3d);

        // calculate new residue
        residue_new = calPointResidue(tr3d._pt_center, region_list, imgAug_list, otf);

        // calculate grad
        for (int i = 0; i < 3; i ++)
        {
            grad[i] = (residue_new - residue_old) / shift3d[i];
        }

        // update shift 3d
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
        std::cout << "step = " << record_step << ", shift3d.mag = " << shift3d.norm() << ", residue_orig = " << residue_orig << ", residue_new = " << residue_new <<  std::endl;
    }
    else
    {
        tr3d._pt_center = pt_center_orig;
        residue_new = residue_orig;
    }
    
    return residue_new;
}


void Shake::updateImgAugList (std::vector<Image>& imgAug_list, std::vector<PixelRange>& region_list, Tracer3D const& tr3d)
{
    // update Augemented image and search range
    for (int id = 0; id < _n_cam_use; id ++)
    {
        // Get the new range based on the updated 3d position
        // TODO: what if the new projected 2d pt is out of image? (n_col_new==1 or n_row_new==1)
        PixelRange region = findRegion(id, tr3d._tr2d_list[id]._pt_center[1], tr3d._tr2d_list[id]._pt_center[0], tr3d._tr2d_list[id]._r_px);
        int n_col_new = region.getNumOfCol();
        int n_row_new = region.getNumOfRow();
        int n_col_old = region_list[id].getNumOfCol();
        int n_row_old = region_list[id].getNumOfRow();

        // save the orig imgAug and update imgAug based on new range
        Image imgAug_old = imgAug_list[id];
        if (n_col_new != n_col_old || 
            n_row_new != n_row_old)
        {
            imgAug_list[id] = Image(n_row_new, n_col_new, 0);
        }
        
        // record the old range
        int row_min_old = region_list[id].row_min;
        int row_max_old = region_list[id].row_max;
        int col_min_old = region_list[id].col_min;
        int col_max_old = region_list[id].col_max;

        // update search range
        region_list[id] = region;

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
                    imgAug_list[id](i,j) = imgAug_old(row-row_min_old, col-col_min_old);
                }
                else
                {
                    imgAug_list[id](i,j) = std::max(0.0, _imgRes_list[id](row, col));
                }
                j ++;
            }
            i ++;
        }   
    }    
}


// Try with image cross correlation
double Shake::calPointResidue (Pt3D const& pt3d, std::vector<PixelRange> const& region_list, std::vector<Image> const& imgAug_list, OTF const& otf)
{
    double residue = 0;
    int cam_id;
    int npts = 0;
    Pt2D pt2d;

    for (int id = 0; id < _n_cam_use; id ++)
    {
        cam_id = _cam_list.useid_list[id];
        pt2d = _cam_list.cam_list[cam_id].project(pt3d);

        // get otf param
        std::vector<double> otf_param = otf.getOTFParam(cam_id, pt3d);

        // Calculating the residual for updated 3D position
        //  residue = sum(sum((I_p - I_o)^2) / npts, all cam) / n_cam
        int i = 0; 
        for (int row = region_list[id].row_min; row < region_list[id].row_max; row ++)
        {
            int j = 0; 
            for (int col = region_list[id].col_min; col < region_list[id].col_max; col ++)
            {
                residue += std::pow(
                    imgAug_list[id](i,j) - 
                    gaussIntensity(col, row, pt2d, otf_param),
                    2
                );

                j ++;
                npts ++;
            }
            i ++;
        }
    }

    residue /= npts;
    return residue;
}


double Shake::calTracerScore (Tracer3D const& tr3d, std::vector<PixelRange> const& region_list, std::vector<Image> const& imgAug_list, OTF const& otf, double score)
{
    // update particle intensity 
    // (ignoring the camera with highest local peak intensity)
    // Calculate intensity for each camera: 
    //  if the tracer is out of the image, directly set the intensity as 0
    //  then calculate the new score without the outlier
    
    std::vector<double> numerator_list(_n_cam_use, 0);
    std::vector<double> denominator_list(_n_cam_use, 0);
    int n_outRange = 0;
    for (int id = 0; id < _n_cam_use; id ++)
    {
        std::vector<double> otf_param = otf.getOTFParam(
            _cam_list.useid_list[id],
            tr3d._pt_center
        );

        // if the tracer is out of the image, directly set the intensity as 0
        if (region_list[id].getNumOfCol() == 1 || 
            region_list[id].getNumOfRow() == 1)
        {
            n_outRange ++;
            if (_n_cam_use-n_outRange < 2)
            {
                return 0;
            } 
            // return 0;
        }

        int i = 0;
        for (int row = region_list[id].row_min; row < region_list[id].row_max; row ++)
        {
            int j = 0;
            for (int col = region_list[id].col_min; col < region_list[id].col_max; col ++)
            {
                numerator_list[id] += imgAug_list[id](i, j); 
                
                denominator_list[id] += gaussIntensity(col, row, tr3d._tr2d_list[id]._pt_center, otf_param);

                j ++;
            }
            i ++;
        }
    }


    // calculate numerator and denominator
    double numerator = 0.0;     
    double denominator = 0.0;
    double threshold = 50; // TODO: add input for this
    int n_notfind = 0;
    for (int id = 0; id < _n_cam_use; id ++)
    {
        if (numerator_list[id] < threshold)
        {
            n_notfind ++;
            continue;
        }
        numerator += numerator_list[id];
        denominator += denominator_list[id];   
    }

    if (_n_cam_use - n_notfind < 2)
    {
        return 0;
    }

    if (denominator < SMALLNUMBER)
    {
        // // for debug
        // std::cerr << "Warning: Shake::calTracerIntensity warning at line " << __LINE__ << ":\n" 
        //           << "denominator is 0!"
        //           << std::endl;
        denominator = SMALLNUMBER;
    }

    // double score_new = score * std::fabs(numerator/denominator);
    double score_new = std::fabs(numerator/denominator);
    return score_new;
}


void Shake::removeGhost(std::vector<Tracer3D>& tr3d_list)
{
    int n_tr3d = tr3d_list.size();

    // compute mean score without zeros
    double mean = 0; 
    int n_mean = 0;
    for (int i = 0; i < n_tr3d; i ++)
    {
        if (_score_list[i] > SMALLNUMBER)
        {
            mean += _score_list[i];
            n_mean ++;
        }
    }
    mean /= n_mean;

    // std::cout << "mean score: " << mean << std::endl;

    // remove ghost tracers if the score is less than _min_score*mean
    // but _score_list is erased
    std::vector<Tracer3D> tr3d_list_new;
    for (int i = n_tr3d-1; i > -1; i --)
    {
        if (_score_list[i] < _score_min * mean)
        {
            tr3d_list.erase(tr3d_list.begin() + i);
            _objID_remove.push_back(i);
            _n_ghost ++;
        }
        else
        {
            _objID_keep.push_back(i);
        }
    }

    // reverse _objID_keep back
    std::reverse(_objID_keep.begin(), _objID_keep.end());
    std::reverse(_objID_remove.begin(), _objID_remove.end());


    std::cout << "\tShake::removeGhost: " << _n_ghost << " ghost tracers are removed." << std::endl;
}


void Shake::removeGhostResidue(std::vector<Tracer3D>& tr3d_list)
{
    int n_tr3d = tr3d_list.size();

    // compute mean score without zeros
    double mean = 0; 
    int n_mean = 0;
    for (int i = 0; i < n_tr3d; i ++)
    {
        mean += _score_list[i];
        n_mean ++;
    }
    mean /= n_mean;

    double std = 0;
    for (int i = 0; i < n_tr3d; i ++)
    {
        std += std::pow(_score_list[i] - mean, 2);
    }
    std = std::sqrt(std / n_mean);

    std::cout << "mean score: " << mean << std::endl;
    std::cout << "std score: " << std << std::endl;

    // remove ghost tracers if the score is less than _min_score*mean
    // but _score_list is erased
    std::vector<Tracer3D> tr3d_list_new;
    for (int i = n_tr3d-1; i > -1; i --)
    {
        // if (_score_list[i] > mean + _score_min * std)
        if (_score_list[i] > _score_min)
        {
            tr3d_list.erase(tr3d_list.begin() + i);
            _objID_remove.push_back(i);
            _n_ghost ++;
        }
        else
        {
            _objID_keep.push_back(i);
        }
    }

    // reverse _objID_keep back to original order: small to large
    std::reverse(_objID_keep.begin(), _objID_keep.end());
    std::reverse(_objID_remove.begin(), _objID_remove.end());

    std::cout << "\tShake::removeGhost: " << _n_ghost << " ghost tracers are removed." << std::endl;
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