#include "Camera.h"

Camera::Camera () {};

Camera::Camera (const Camera& c)
    : _type(c._type), _pinhole_param(c._pinhole_param), _poly_param(c._poly_param) {}

Camera::Camera(std::istream& is)
{
    loadParameters(is);
}

Camera::Camera(std::string file_name)
{
    loadParameters(file_name);
}

void Camera::loadParameters (std::istream& is)
{
    std::string type_name;
    is >> type_name;
    if (type_name == "PINHOLE")
    {
        _type = PINHOLE;

        // do not read errors
        std::string useless;
        is >> useless;
        is >> useless;
        
        // read camera matrix
        _pinhole_param.cam_mtx = Matrix<double>(3,3,is);

        // initialize distortion parameters
        _pinhole_param.dist_coeff.clear();
        _pinhole_param.is_distorted = false;
        _pinhole_param.n_dist_coeff = 0;

        // read distortion coefficients
        std::string dist_coeff_str;
        is >> dist_coeff_str;
        std::stringstream dist_coeff_stream(dist_coeff_str);
        double dist_coeff;
        int id = 0;
        while (dist_coeff_stream >> dist_coeff)
        {
            _pinhole_param.dist_coeff.push_back(dist_coeff);
            id ++;

            if (dist_coeff > SMALLNUMBER)
            {
                _pinhole_param.is_distorted = true;
            }

            if (dist_coeff_stream.peek() == ',')
            {
                dist_coeff_stream.ignore();
            } 
        }
        if (_pinhole_param.is_distorted)
        {
            _pinhole_param.n_dist_coeff = id;
            if (id != 4 && id != 5 && id != 8 && id != 12)
            {
                // Note: id == 14
                // additional distortion by projecting onto a tilt plane
                // not supported yet

                std::cerr << "Camera::LoadParameters line " << __LINE__ << " : Error: number of distortion coefficients is wrong: " << id << std::endl;
                throw error_type;
            }
        }

        // do not read rotation vector
        is >> useless;

        // read rotation matrix
        _pinhole_param.r_mtx = Matrix<double>(3,3,is);

        // read inverse rotation matrix
        _pinhole_param.r_mtx_inv = Matrix<double>(3,3,is);

        // read translation vector
        _pinhole_param.t_vec = Pt3D(is);

        // read inverse translation vector
        _pinhole_param.t_vec_inv = Pt3D(is);
    }
    else if (type_name == "POLYNOMIAL")
    {
        _type = POLYNOMIAL;
        
        // do not read errors
        std::string useless;
        is >> useless;

        // read reference plane
        std::string ref_plane_str;
        is >> ref_plane_str;
        std::stringstream plane_stream(ref_plane_str);
        std::string temp;
        std::getline(plane_stream, temp, ',');
        int derivative_id[2] = {1,2};
        if (temp == "REF_X")
        {
            _poly_param.ref_plane = REF_X;
            derivative_id[0] = 2;
            derivative_id[1] = 3;
        }
        else if (temp == "REF_Y")
        {
            _poly_param.ref_plane = REF_Y;
            derivative_id[0] = 1;
            derivative_id[1] = 3;
        }
        else if (temp == "REF_Z")
        {
            _poly_param.ref_plane = REF_Z;
            derivative_id[0] = 1;
            derivative_id[1] = 2;
        }
        else
        {
            std::cerr << "Camera::LoadParameters line " << __LINE__ << " : Error: reference plane is wrong: " << temp << std::endl;
            throw error_type;
        }

        std::getline(plane_stream, temp, ',');
        _poly_param.plane[0] = std::stod(temp);
        std::getline(plane_stream, temp, ',');
        _poly_param.plane[1] = std::stod(temp); 

        // read number of coefficients
        is >> _poly_param.n_coeff;

        // read u coefficients
        _poly_param.u_coeffs = Matrix<double>(_poly_param.n_coeff,4,is);

        // read v coefficients
        _poly_param.v_coeffs = Matrix<double>(_poly_param.n_coeff,4,is);

        // calculate du, dv coefficients
        _poly_param.du_coeffs = Matrix<double>(_poly_param.n_coeff*2,4,0);
        _poly_param.dv_coeffs = Matrix<double>(_poly_param.n_coeff*2,4,0);
        for (int i = 0; i < _poly_param.n_coeff; i ++)
        {
            // calculate du coeff
            _poly_param.du_coeffs(i,0) = _poly_param.u_coeffs(i,0) * _poly_param.u_coeffs(i,derivative_id[0]);   
            _poly_param.du_coeffs(i,_poly_param.ref_plane) = _poly_param.u_coeffs(i,_poly_param.ref_plane);
            _poly_param.du_coeffs(i,derivative_id[0]) = std::max(_poly_param.u_coeffs(i,derivative_id[0]) - 1, 0.0);
            _poly_param.du_coeffs(i,derivative_id[1]) = _poly_param.u_coeffs(i,derivative_id[1]);

            int j = i + _poly_param.n_coeff;
            _poly_param.du_coeffs(j,0) = _poly_param.u_coeffs(i,0) * _poly_param.u_coeffs(i,derivative_id[1]);
            _poly_param.du_coeffs(j,_poly_param.ref_plane) = _poly_param.u_coeffs(i,_poly_param.ref_plane);
            _poly_param.du_coeffs(j,derivative_id[0]) = _poly_param.u_coeffs(i,derivative_id[0]);
            _poly_param.du_coeffs(j,derivative_id[1]) = std::max(_poly_param.u_coeffs(i,derivative_id[1]) - 1, 0.0);

            // calculate dv coeff
            _poly_param.dv_coeffs(i,0) = _poly_param.v_coeffs(i,0) * _poly_param.v_coeffs(i,derivative_id[0]);
            _poly_param.dv_coeffs(i,_poly_param.ref_plane) = _poly_param.v_coeffs(i,_poly_param.ref_plane);
            _poly_param.dv_coeffs(i,derivative_id[0]) = std::max(_poly_param.v_coeffs(i,derivative_id[0]) - 1, 0.0);
            _poly_param.dv_coeffs(i,derivative_id[1]) = _poly_param.v_coeffs(i,derivative_id[1]);

            _poly_param.dv_coeffs(j,0) = _poly_param.v_coeffs(i,0) * _poly_param.v_coeffs(i,derivative_id[1]);
            _poly_param.dv_coeffs(j,_poly_param.ref_plane) = _poly_param.v_coeffs(i,_poly_param.ref_plane);
            _poly_param.dv_coeffs(j,derivative_id[0]) = _poly_param.v_coeffs(i,derivative_id[0]);
            _poly_param.dv_coeffs(j,derivative_id[1]) = std::max(_poly_param.v_coeffs(i,derivative_id[1]) - 1, 0.0);
        }
    }
    else 
    {
        std::cerr << "Camera::LoadParameters line " << __LINE__ << " : Error: unknown camera type: " << type_name << std::endl;
        throw error_type;
    }
}

void Camera::loadParameters (std::string file_name)
{
    std::ifstream infile(file_name.c_str(), std::ios::in);
    std::string line;
    std::stringstream file_content;
    while (std::getline(infile, line)) {
        size_t comment_pos = line.find('#');
        if (comment_pos > 0) 
        {
            if (comment_pos < std::string::npos) 
            {
                line.erase(comment_pos);
    		}
    	}
        else if (comment_pos == 0)
        {
            continue;
        }

        file_content << line << '\t';
    }
    infile.close();

    loadParameters(file_content);
}

Pt3D Camera::rmtxTorvec (Matrix<double> const& r_mtx)
{
    Pt3D r_vec;

    double tr = (myMATH::trace<double>(r_mtx) - 1) / 2;
    tr = tr > 1 ? 1 : tr < -1 ? -1 : tr;
    double theta = std::acos(tr);
    double s = std::sin(theta);

    if (s > SMALLNUMBER)
    {
        double ratio = theta / (2 * s);
        r_vec[0] = (r_mtx(2,1) - r_mtx(1,2)) * ratio;
        r_vec[1] = (r_mtx(0,2) - r_mtx(2,0)) * ratio;
        r_vec[2] = (r_mtx(1,0) - r_mtx(0,1)) * ratio;
    }
    else if (tr > 0)
    {
        r_vec[0] = 0;
        r_vec[1] = 0;
        r_vec[2] = 0;
    }
    else
    {
        r_vec[0] = theta * std::sqrt((r_mtx(0,0) + 1) / 2);
        r_vec[1] = theta * std::sqrt((r_mtx(1,1) + 1) / 2) * (r_mtx(0,1) > 0 ? 1 : -1);
        r_vec[2] = theta * std::sqrt((r_mtx(2,2) + 1) / 2) * (r_mtx(0,2) > 0 ? 1 : -1);
    }

    return r_vec;
}

void Camera::saveParameters (std::string file_name)
{
    std::ofstream outfile(file_name.c_str(), std::ios::out);
    if (_type == PINHOLE)
    {
        outfile << "# Camera Model: (PINHOLE/POLYNOMIAL)" << std::endl;
        outfile << "PINHOLE" << std::endl;
        outfile << "# Camera Calibration Error: \nNone\n# Pose Calibration Error: \nNone" << std::endl;
        
        outfile << "# Camera Matrix: " << std::endl;
        _pinhole_param.cam_mtx.write(outfile);
        
        outfile << "# Distortion Coefficients: " << std::endl;
        int size = _pinhole_param.dist_coeff.size();
        for (int i = 0; i < size-1; i ++)
        {
            outfile << _pinhole_param.dist_coeff[i] << ",";
        }
        outfile << _pinhole_param.dist_coeff[size-1] << std::endl;

        outfile << "# Rotation Vector: " << std::endl;
        Pt3D r_vec = rmtxTorvec(_pinhole_param.r_mtx);
        r_vec.transpose().write(outfile);

        outfile << "# Rotation Matrix: " << std::endl;
        _pinhole_param.r_mtx.write(outfile);

        outfile << "# Inverse of Rotation Matrix: " << std::endl;
        _pinhole_param.r_mtx_inv.write(outfile);

        outfile << "# Translation Vector: " << std::endl;
        _pinhole_param.t_vec.transpose().write(outfile);

        outfile << "# Inverse of Translation Vector: " << std::endl;
        _pinhole_param.t_vec_inv.transpose().write(outfile);

        outfile.close();
    }
    else if (_type == POLYNOMIAL)
    {
        outfile << "# Camera Model: (PINHOLE/POLYNOMIAL)" << std::endl;
        outfile << "POLYNOMIAL" << std::endl;
        outfile << "# Camera Calibration Error: \nNone" << std::endl;

        outfile << "# Reference Plane: (REF_X/REF_Y/REF_Z,coordinate,coordinate)" << std::endl;
        if (_poly_param.ref_plane == REF_X)
        {
            outfile << "REF_X,";
        }
        else if (_poly_param.ref_plane == REF_Y)
        {
            outfile << "REF_Y,";
        }
        else if (_poly_param.ref_plane == REF_Z)
        {
            outfile << "REF_Z,";
        }
        outfile << _poly_param.plane[0] << "," << _poly_param.plane[1] << std::endl;

        outfile << "# Number of Coefficients: " << std::endl;
        outfile << _poly_param.n_coeff << std::endl;

        outfile << "# U_Coeff,X_Power,Y_Power,Z_Power" << std::endl;
        _poly_param.u_coeffs.write(outfile);

        outfile << "# V_Coeff,X_Power,Y_Power,Z_Power" << std::endl;
        _poly_param.v_coeffs.write(outfile);

        outfile.close();
    }
    else
    {
        std::cerr << "Camera::SaveParameters line " << __LINE__ << " : Error: unknown camera type: " << _type << std::endl;
        throw error_type;
    }
}

//                                                 
// Project world coordinate [mm] to image points [px] 
//
Pt2D Camera::project (Pt3D const& pt_world)
{
    if (_type == PINHOLE)
    {
        return distort(worldToUndistImg(pt_world));
    }
    else if (_type == POLYNOMIAL)
    {
        return polyProject(pt_world);
    }
    else
    {
        std::cerr << "Camera::Project line " << __LINE__ << " : Error: unknown camera type: " << _type << std::endl;
        throw error_type;
    }
}

// Pinhole  model
Pt2D Camera::worldToUndistImg (Pt3D const& pt_world)
{
    Pt3D temp = _pinhole_param.r_mtx * pt_world + _pinhole_param.t_vec;

    Pt2D pt_img_mm;
    pt_img_mm[0] = temp[2] ? (temp[0]/temp[2]) : temp[0];
    pt_img_mm[1] = temp[2] ? (temp[1]/temp[2]) : temp[1];

    return pt_img_mm;
}

Pt2D Camera::distort (Pt2D const& pt_img_undist)
{
    // opencv distortion model
    double x = pt_img_undist[0];
    double y = pt_img_undist[1];
    double xd, yd;

    // judge if distortion is needed
    if (_pinhole_param.is_distorted)
    {
        double r2 = x*x + y*y;
        double r4 = r2*r2;
        double r6 = r4*r2;

        double a1 = 2 * x*y;
        double a2 = r2 + 2 * x*x;
        double a3 = r2 + 2 * y*y;

        double cdist = 1 + _pinhole_param.dist_coeff[0]*r2 + _pinhole_param.dist_coeff[1]*r4;
        if (_pinhole_param.n_dist_coeff > 4)
        {
            cdist += _pinhole_param.dist_coeff[4]*r6;
        }         

        double icdist2 = 1.0;
        if (_pinhole_param.n_dist_coeff > 5)
        {
            icdist2 = 1.0 / (1 + _pinhole_param.dist_coeff[5]*r2 + _pinhole_param.dist_coeff[6]*r4 + _pinhole_param.dist_coeff[7]*r6);
        }
        // double icdist2 = 1.0 / (1 + _pinhole_param.dist_coeff[5]*r2 + _pinhole_param.dist_coeff[6]*r4 + _pinhole_param.dist_coeff[7]*r6);

        xd = x*cdist*icdist2 + _pinhole_param.dist_coeff[2]*a1 + _pinhole_param.dist_coeff[3]*a2;
        yd = y*cdist*icdist2 + _pinhole_param.dist_coeff[2]*a3 + _pinhole_param.dist_coeff[3]*a1;
        if (_pinhole_param.n_dist_coeff > 8)
        {
            // additional distortion by projecting onto a tilt plane (14)
            // not supported yet
            xd += _pinhole_param.dist_coeff[8] *r2 + _pinhole_param.dist_coeff[9] *r4;
            yd += _pinhole_param.dist_coeff[10]*r2 + _pinhole_param.dist_coeff[11]*r4;
        }
    }
    else
    {
        xd = x;
        yd = y;
    }
    
    Pt2D pt_img_dist(
        xd * _pinhole_param.cam_mtx(0,0) + _pinhole_param.cam_mtx(0,2),
        yd * _pinhole_param.cam_mtx(1,1) + _pinhole_param.cam_mtx(1,2)
    );

    return pt_img_dist;
}

// Polynomial model
Pt2D Camera::polyProject (Pt3D const& pt_world)
{
    double u = 0;
    double v = 0;

    for (int i = 0; i < _poly_param.u_coeffs.getDimRow(); i ++)
    {
        double u_val = _poly_param.u_coeffs(i,0);
        double v_val = _poly_param.v_coeffs(i,0);

        for (int j = 1; j < 4; j ++)
        {
            u_val *= std::pow(pt_world[j-1], _poly_param.u_coeffs(i,j));
            v_val *= std::pow(pt_world[j-1], _poly_param.v_coeffs(i,j));
        }

        u += u_val;
        v += v_val;
    }

    return Pt2D(u,v);
}


//
// Calculate line of sight
//
Line3D Camera::lineOfSight (Pt2D const& pt_img_dist)
{
    if (_type == PINHOLE)
    {
        return pinholeLine(undistort(pt_img_dist));
    }
    else if (_type == POLYNOMIAL)
    {
        return polyLineOfSight(pt_img_dist);
    }
    else
    {
        std::cerr << "Camera::LineOfSight line " << __LINE__ << " : Error: unknown camera type: " << _type << std::endl;
        throw error_type;
    }
}

// Pinhole  model
Pt2D Camera::undistort (Pt2D const& pt_img_dist)
{
    double fx = _pinhole_param.cam_mtx(0,0);
    double fy = _pinhole_param.cam_mtx(1,1);
    double cx = _pinhole_param.cam_mtx(0,2);
    double cy = _pinhole_param.cam_mtx(1,2);

    double u = pt_img_dist[0];
    double v = pt_img_dist[1];
    double x = (u - cx) / fx;
    double y = (v - cy) / fy;
    Pt2D img_undist(x,y);
    Pt2D img_dist(0,0);

    if (_pinhole_param.is_distorted)
    {
        double x0 = x;
        double y0 = y;
        double error = std::numeric_limits<double>::max();

        // undistort iteratively
        int iter = 0;
        while (error > UNDISTORT_EPS && iter < UNDISTORT_MAX_ITER)
        {   
            double r2 = x*x + y*y;
            double a1 = 2 * x*y;
            double a2 = r2 + 2 * x*x;
            double a3 = r2 + 2 * y*y;

            // double icdist = (1 + ( (k[7]*r2+k[6])*r2 + k[5] )*r2)
            //                /(1 + ( (k[4]*r2+k[1])*r2 + k[0] )*r2);
            double icdist;
            if (_pinhole_param.n_dist_coeff == 4)
            {
                icdist = 1 / (1 + (_pinhole_param.dist_coeff[1]*r2 + _pinhole_param.dist_coeff[0])*r2);
            }
            else if (_pinhole_param.n_dist_coeff == 5)
            {
                icdist = 1 / (1 + ((_pinhole_param.dist_coeff[4]*r2 + _pinhole_param.dist_coeff[1])*r2 + _pinhole_param.dist_coeff[0])*r2);
            }
            else
            {
                icdist = (1 + ((_pinhole_param.dist_coeff[7]*r2 + _pinhole_param.dist_coeff[6])*r2 + _pinhole_param.dist_coeff[5])*r2)/(1 + ((_pinhole_param.dist_coeff[4]*r2 + _pinhole_param.dist_coeff[1])*r2 + _pinhole_param.dist_coeff[0])*r2);
            }

            if (icdist < 0)
            {
                img_undist[0] = (u - cx) / fx;
                img_undist[1] = (v - cy) / fy;
                break;
            }

            // update x, y
            double deltaX = _pinhole_param.dist_coeff[2] * a1 + _pinhole_param.dist_coeff[3] * a2;
            double deltaY = _pinhole_param.dist_coeff[2] * a3 + _pinhole_param.dist_coeff[3] * a1;
            if (_pinhole_param.n_dist_coeff > 8)
            {
                deltaX += _pinhole_param.dist_coeff[8] * r2 + _pinhole_param.dist_coeff[9] * r2 * r2;
                deltaY += _pinhole_param.dist_coeff[10] * r2 + _pinhole_param.dist_coeff[11] * r2 * r2;
            }
            x = (x0 - deltaX) * icdist;
            y = (y0 - deltaY) * icdist;

            // calculate distorted xd, yd
            img_undist[0] = x;
            img_undist[1] = y;
            img_dist = distort(img_undist);

            // update error
            error = std::sqrt(std::pow(img_dist[0]-u, 2) + std::pow(img_dist[1]-v, 2));

            iter ++;
        }
    }

    return img_undist;
}

Line3D Camera::pinholeLine (Pt2D const& pt_img_undist)
{
    Pt3D pt_world(pt_img_undist[0], pt_img_undist[1], 1.0);

    // calculate line of sight
    pt_world = _pinhole_param.r_mtx_inv * pt_world + _pinhole_param.t_vec_inv;
    
    Pt3D unit_vec = myMATH::createUnitVector(_pinhole_param.t_vec_inv, pt_world);

    Line3D line = {_pinhole_param.t_vec_inv, unit_vec};

    return line;
}

// Polynomial model
Pt3D Camera::polyImgToWorld (Pt2D const& pt_img_dist, double plane_world)
{
    Pt3D pt_world(
        (double)rand() / RAND_MAX,
        (double)rand() / RAND_MAX,
        (double)rand() / RAND_MAX
    );

    switch (_poly_param.ref_plane)
    {
        case REF_X:
            pt_world[0] = plane_world;
            break;
        case REF_Y:
            pt_world[1] = plane_world;
            break;
        case REF_Z:
            pt_world[2] = plane_world;
            break;
        default:
            std::cerr << "Camera::PolyImgToWorld line " << __LINE__ << " : Error: unknown reference plane: " << _poly_param.ref_plane << std::endl;
            throw error_type;
    }

    Matrix<double> jacobian(2,2,0);
    Matrix<double> jacobian_inv(2,2,0);
    Pt2D pt_img_temp;
    double du,dv,dx,dy;
    double err = std::numeric_limits<double>::max();
    int iter = 0;
    jacobian = Matrix<double>(2,2,0);
    jacobian_inv = Matrix<double>(2,2,0);

    while (err > UNDISTORT_EPS && iter < UNDISTORT_MAX_ITER)
    {
        // calculate jacobian matrix (e.g. REF_Z)
        // du = u_true - u, dv = v_true - v
        // |du| = | du/dx, du/dy | |dx|
        // |dv| = | dv/dx, dv/dy | |dy|
        // x = x0 + dx 
        // y = y0 + dy
        jacobian *= 0;
        jacobian_inv *= 0;
        pt_img_temp = polyProject(pt_world);
        du = pt_img_dist[0] - pt_img_temp[0];
        dv = pt_img_dist[1] - pt_img_temp[1];

        // calculate jacobian matrix    
        for (int i = 0; i < _poly_param.n_coeff; i ++)
        {
            double dudx = _poly_param.du_coeffs(i,0);
            double dudy = _poly_param.du_coeffs(i+_poly_param.n_coeff,0);
            double dvdx = _poly_param.dv_coeffs(i,0);
            double dvdy = _poly_param.dv_coeffs(i+_poly_param.n_coeff,0);

            for (int j = 1; j < 4; j ++)
            {
                dudx *= std::pow(pt_world[j-1], _poly_param.du_coeffs(i,j));
                dudy *= std::pow(pt_world[j-1], _poly_param.du_coeffs(i+_poly_param.n_coeff,j));

                dvdx *= std::pow(pt_world[j-1], _poly_param.dv_coeffs(i,j));
                dvdy *= std::pow(pt_world[j-1], _poly_param.dv_coeffs(i+_poly_param.n_coeff,j));
            }

            jacobian(0,0) += dudx;
            jacobian(0,1) += dudy;
            jacobian(1,0) += dvdx;
            jacobian(1,1) += dvdy;
        }

        // calculate dx, dy
        jacobian_inv = myMATH::inverse(jacobian, "det");
        pt_img_temp = polyProject(pt_world);
        du = pt_img_dist[0] - pt_img_temp[0];
        dv = pt_img_dist[1] - pt_img_temp[1];

        dx = jacobian_inv(0,0) * du + jacobian_inv(0,1) * dv;
        dy = jacobian_inv(1,0) * du + jacobian_inv(1,1) * dv;

        // update pt_world
        switch (_poly_param.ref_plane)
        {
            case REF_X:
                pt_world[1] += dx;
                pt_world[2] += dy;
                break;
            case REF_Y:
                pt_world[0] += dx;
                pt_world[2] += dy;
                break;
            case REF_Z:
                pt_world[0] += dx;
                pt_world[1] += dy;
                break;
            default:
                std::cerr << "Camera::PolyImgToWorld line " << __LINE__ << " : Error: unknown reference plane: " << _poly_param.ref_plane << std::endl;
                throw error_type;
        }

        // update error, iter
        err = std::sqrt(du*du + dv*dv);
        iter ++;
    }

    // for debug
    // std::cout << "iter: " << iter << std::endl;
    // std::cout << "err: " << err << std::endl;

    return pt_world;
}

Line3D Camera::polyLineOfSight (Pt2D const& pt_img_dist)
{
    Pt3D pt_world_1 = polyImgToWorld(pt_img_dist, _poly_param.plane[0]);
    Pt3D pt_world_2 = polyImgToWorld(pt_img_dist, _poly_param.plane[1]);
    Pt3D unit_vec = myMATH::createUnitVector(pt_world_1, pt_world_2);

    Line3D line = {pt_world_1, unit_vec};

    return line;
}