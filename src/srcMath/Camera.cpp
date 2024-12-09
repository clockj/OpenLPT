#include "Camera.h"

Camera::Camera () {};

Camera::Camera (const Camera& c)
    : _type(c._type), _pinhole_param(c._pinhole_param), _poly_param(c._poly_param), _pinplate_param(c._pinplate_param) {}

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

        // read image size (n_row,n_col)
        std::string img_size_str;
        is >> img_size_str;
        std::stringstream img_size_stream(img_size_str);
        std::string temp;
        std::getline(img_size_stream, temp, ',');
        _pinhole_param.n_row = std::stoi(temp);
        std::getline(img_size_stream, temp, ',');
        _pinhole_param.n_col = std::stoi(temp);
        
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

        // read image size (n_row,n_col)
        std::string img_size_str;
        is >> img_size_str;
        std::stringstream img_size_stream(img_size_str);
        std::string temp;
        std::getline(img_size_stream, temp, ',');
        _poly_param.n_row = std::stoi(temp);
        std::getline(img_size_stream, temp, ',');
        _poly_param.n_col = std::stoi(temp);

        // read reference plane
        std::string ref_plane_str;
        is >> ref_plane_str;
        std::stringstream plane_stream(ref_plane_str);
        // std::string temp;
        std::getline(plane_stream, temp, ',');
        if (temp == "REF_X")
        {
            _poly_param.ref_plane = REF_X;
        }
        else if (temp == "REF_Y")
        {
            _poly_param.ref_plane = REF_Y;
        }
        else if (temp == "REF_Z")
        {
            _poly_param.ref_plane = REF_Z;
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
        updatePolyDuDv ();
    }
    else if (type_name == "PINPLATE")
    {
        _type = PINPLATE;

        // do not read errors
        std::string useless;
        is >> useless;
        is >> useless;

        // read image size (n_row,n_col)
        std::string img_size_str;
        is >> img_size_str;
        std::stringstream img_size_stream(img_size_str);
        std::string temp;
        std::getline(img_size_stream, temp, ',');
        _pinplate_param.n_row = std::stoi(temp);
        std::getline(img_size_stream, temp, ',');
        _pinplate_param.n_col = std::stoi(temp);
        
        // read camera matrix
        _pinplate_param.cam_mtx = Matrix<double>(3,3,is);

        // initialize distortion parameters
        _pinplate_param.dist_coeff.clear();
        _pinplate_param.is_distorted = false;
        _pinplate_param.n_dist_coeff = 0;

        // read distortion coefficients
        std::string dist_coeff_str;
        is >> dist_coeff_str;
        std::stringstream dist_coeff_stream(dist_coeff_str);
        double dist_coeff;
        int id = 0;
        while (dist_coeff_stream >> dist_coeff)
        {
            _pinplate_param.dist_coeff.push_back(dist_coeff);
            id ++;

            if (dist_coeff > SMALLNUMBER)
            {
                _pinplate_param.is_distorted = true;
            }

            if (dist_coeff_stream.peek() == ',')
            {
                dist_coeff_stream.ignore();
            } 
        }
        if (_pinplate_param.is_distorted)
        {
            _pinplate_param.n_dist_coeff = id;
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
        _pinplate_param.r_mtx = Matrix<double>(3,3,is);

        // read inverse rotation matrix
        _pinplate_param.r_mtx_inv = Matrix<double>(3,3,is);

        // read translation vector
        _pinplate_param.t_vec = Pt3D(is);

        // read inverse translation vector
        _pinplate_param.t_vec_inv = Pt3D(is);

        // read reference point of refractive plate
        _pinplate_param.plane.pt = Pt3D(is);

        // read normal vector of refractive plate
        _pinplate_param.plane.norm_vector = Pt3D(is);

        // read refractive index array
        std::string refract_str;
        is >> refract_str;
        std::stringstream refract_stream(refract_str);
        double refract;
        while (refract_stream >> refract)
        {
            _pinplate_param.refract_array.push_back(refract);
            if (refract_stream.peek() == ',')
            {
                refract_stream.ignore();
            } 
        }

        // read width of the refractive plate
        std::string w_str;
        is >> w_str;
        std::stringstream w_stream(w_str);
        double w;
        while (w_stream >> w)
        {
            _pinplate_param.w_array.push_back(w);
            if (w_stream.peek() == ',')
            {
                w_stream.ignore();
            } 
        }

        // check number of plate is consistent
        if (_pinplate_param.refract_array.size() != _pinplate_param.w_array.size()+2 || _pinplate_param.refract_array.size() < 3)
        {
            std::cerr << "Camera::LoadParameters line " << __LINE__ 
                      << " : Error: number of refractive index and width is not consistent (n_plate_min=1): " 
                      << _pinplate_param.refract_array.size() << " vs " << _pinplate_param.w_array.size() << std::endl;
            throw error_type;
        }
        _pinplate_param.n_plate = _pinplate_param.w_array.size();

        // update closest point to the camera
        _pinplate_param.pt3d_closest = _pinplate_param.plane.pt;
        for (int i = 0; i < _pinplate_param.n_plate; i ++)
        {
            _pinplate_param.pt3d_closest = _pinplate_param.pt3d_closest - _pinplate_param.plane.norm_vector * _pinplate_param.w_array[i];
        }

        // read projection tolerance
        double tol;
        is >> tol;
        _pinplate_param.proj_tol2 = tol * tol;

        // read maximum number of iterations for projection
        is >> _pinplate_param.proj_nmax;

        // read learning rate if exists
        if (is >> _pinplate_param.lr)
        {
            // do nothing
        }
        else
        {
            _pinplate_param.lr = 0.1;
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

    if (!infile)
    {
        std::cerr << "Camera::LoadParameters line " << __LINE__ << " : Error: cannot open file: " << file_name << std::endl;
        throw error_io;
    }

    std::string line;
    std::stringstream file_content;
    while (std::getline(infile, line)) 
    {
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

void Camera::updatePolyDuDv ()
{
    if (_type == POLYNOMIAL)
    {
        int derivative_id[2] = {1,2};
        if (_poly_param.ref_plane == REF_X)
        {
            derivative_id[0] = 2;
            derivative_id[1] = 3;
        }
        else if (_poly_param.ref_plane == REF_Y)
        {
            derivative_id[0] = 1;
            derivative_id[1] = 3;
        }
        else if (_poly_param.ref_plane == REF_Z)
        {
            derivative_id[0] = 1;
            derivative_id[1] = 2;
        }
        else
        {
            std::cerr << "Camera::updatePolyDuDv line " << __LINE__ 
                      << " : Error: reference plane is wrong: " 
                      << _poly_param.ref_plane << std::endl;
            throw error_type;
        }

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

        outfile << "# Image Size: (n_row,n_col)" << std::endl;
        outfile << _pinhole_param.n_row << "," << _pinhole_param.n_col << std::endl;
        
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

        outfile << "# Image Size: (n_row,n_col)" << std::endl;
        outfile << _poly_param.n_row << "," << _poly_param.n_col << std::endl;

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
    else if (_type == PINPLATE)
    {
        outfile << "# Camera Model: (PINHOLE/POLYNOMIAL)" << std::endl;
        outfile << "PINPLATE" << std::endl;
        outfile << "# Camera Calibration Error: \nNone\n# Pose Calibration Error: \nNone" << std::endl;

        outfile << "# Image Size: (n_row,n_col)" << std::endl;
        outfile << _pinplate_param.n_row << "," << _pinplate_param.n_col << std::endl;
        
        outfile << "# Camera Matrix: " << std::endl;
        _pinplate_param.cam_mtx.write(outfile);
        
        outfile << "# Distortion Coefficients: " << std::endl;
        int size = _pinplate_param.dist_coeff.size();
        for (int i = 0; i < size-1; i ++)
        {
            outfile << _pinplate_param.dist_coeff[i] << ",";
        }
        outfile << _pinplate_param.dist_coeff[size-1] << std::endl;

        outfile << "# Rotation Vector: " << std::endl;
        Pt3D r_vec = rmtxTorvec(_pinplate_param.r_mtx);
        r_vec.transpose().write(outfile);

        outfile << "# Rotation Matrix: " << std::endl;
        _pinplate_param.r_mtx.write(outfile);

        outfile << "# Inverse of Rotation Matrix: " << std::endl;
        _pinplate_param.r_mtx_inv.write(outfile);

        outfile << "# Translation Vector: " << std::endl;
        _pinplate_param.t_vec.transpose().write(outfile);

        outfile << "# Inverse of Translation Vector: " << std::endl;
        _pinplate_param.t_vec_inv.transpose().write(outfile);

        outfile << "# Reference Point of Refractive Plate: " << std::endl;
        _pinplate_param.plane.pt.transpose().write(outfile);

        outfile << "# Normal Vector of Refractive Plate: " << std::endl;
        _pinplate_param.plane.norm_vector.transpose().write(outfile);

        outfile << "# Refractive Index Array: " << std::endl;
        for (int i = 0; i < _pinplate_param.n_plate+1; i ++)
        {
            outfile << _pinplate_param.refract_array[i] << ",";
        }
        outfile << _pinplate_param.refract_array[_pinplate_param.n_plate+1] << std::endl;

        outfile << "# Width of the Refractive Plate: (same as other physical units)" << std::endl;
        for (int i = 0; i < _pinplate_param.n_plate-1; i ++)
        {
            outfile << _pinplate_param.w_array[i] << ",";
        }
        outfile << _pinplate_param.w_array[_pinplate_param.n_plate-1] << std::endl;

        outfile << "# Projection Tolerance: " << std::endl;
        outfile << std::sqrt(_pinplate_param.proj_tol2) << std::endl;

        outfile << "# Maximum Number of Iterations for Projection: " << std::endl;
        outfile << _pinplate_param.proj_nmax << std::endl;

        outfile << "# Optimization Rate: " << std::endl;
        outfile << _pinplate_param.lr << std::endl;

        outfile.close();
    }
    else
    {
        std::cerr << "Camera::SaveParameters line " << __LINE__ << " : Error: unknown camera type: " << _type << std::endl;
        throw error_type;
    }
}


// 
// Get image size
//
int Camera::getNRow () const
{
    if (_type == PINHOLE)
    {
        return _pinhole_param.n_row;
    }
    else if (_type == POLYNOMIAL)
    {
        return _poly_param.n_row;
    }
    else if (_type == PINPLATE)
    {
        return _pinplate_param.n_row;
    }
    else
    {
        std::cerr << "Camera::GetNRow line " << __LINE__ << " : Error: unknown camera type: " << _type << std::endl;
        throw error_type;
    }
}

int Camera::getNCol () const
{
    if (_type == PINHOLE)
    {
        return _pinhole_param.n_col;
    }
    else if (_type == POLYNOMIAL)
    {
        return _poly_param.n_col;
    }
    else if (_type == PINPLATE)
    {
        return _pinplate_param.n_col;
    }
    else
    {
        std::cerr << "Camera::GetNCol line " << __LINE__ << " : Error: unknown camera type: " << _type << std::endl;
        throw error_type;
    }
}


//                                                 
// Project world coordinate [mm] to image points [px] 
//
Pt2D Camera::project (Pt3D const& pt_world, bool is_print_detail) const
{
    if (_type == PINHOLE)
    {
        return distort(worldToUndistImg(pt_world, _pinhole_param), _pinhole_param);
    }
    else if (_type == POLYNOMIAL)
    {
        return polyProject(pt_world);
    }
    else if (_type == PINPLATE)
    {
        std::tuple<bool,Pt3D,double> result = refractPlate(pt_world);
        if (is_print_detail)
        {
            std::cout << "Camera::Project line " << __LINE__ << " : result (is_parallel,error^2): " << std::get<0>(result) << "," << std::get<2>(result) << std::endl;
        }
        if (std::get<0>(result))
        {
            double value = std::numeric_limits<double>::lowest();
            return Pt2D(value, value);
        }
        // Pt3D pt_refract = std::get<1>(result);
        return distort(worldToUndistImg(std::get<1>(result), _pinplate_param), _pinplate_param);
    }
    else
    {
        std::cerr << "Camera::Project line " << __LINE__ << " : Error: unknown camera type: " << _type << std::endl;
        throw error_type;
    }
}

// Pinhole  model
Pt2D Camera::worldToUndistImg (Pt3D const& pt_world, PinholeParam const& param) const
{
    Pt3D temp = param.r_mtx * pt_world + param.t_vec;

    Pt2D pt_img_mm;
    pt_img_mm[0] = temp[2] ? (temp[0]/temp[2]) : temp[0];
    pt_img_mm[1] = temp[2] ? (temp[1]/temp[2]) : temp[1];

    return pt_img_mm;
}

// Plate refraction model
std::tuple<bool, Pt3D, double> Camera::refractPlate (Pt3D const& pt_world) const
{
    Plane3D plane;
    Line3D line;
    Pt3D pt_cross;
    Pt3D pt_init;
    Pt3D pt_last;
    std::tuple<bool, Pt3D, double> result = {false, {0,0,0}, -1};
    bool is_parallel = false;

    double cos_1; 
    double cos_2;
    double factor;
    double refract_ratio;
    double diff_vec[3] = {0,0,0};
    double delta = -1;
    for (int iter = 1; iter < _pinplate_param.proj_nmax; iter ++)
    {
        plane = _pinplate_param.plane;

        if (iter == 1)
        {
            line.pt = pt_world;
            line.unit_vector = myMATH::createUnitVector(pt_world, _pinplate_param.t_vec_inv);
            is_parallel = myMATH::crossPoint(pt_cross, line, plane);
            if (is_parallel)
            {
                std::get<0>(result) = is_parallel;
                return result;
            }

            pt_init = pt_cross;
        } 
        else
        {
            line.pt = pt_world;
            pt_cross = pt_init;

            line.unit_vector = myMATH::createUnitVector(line.pt, pt_cross);
        }  

        // forward projection: find intersection with each plane
        int plane_id = 0;
        for (plane_id = 1; plane_id < _pinplate_param.n_plate+1; plane_id ++)
        {
            refract_ratio = _pinplate_param.refract_array[plane_id-1] / _pinplate_param.refract_array[plane_id];
            cos_1 = myMATH::dot(line.unit_vector, plane.norm_vector);
            cos_2 = 1 - refract_ratio*refract_ratio*(1 - cos_1*cos_1);
            if (cos_2 <= 0)
            {
                std::cout << "Camera::RefractPlate line " << __LINE__ << " : Error: total internal reflection" << std::endl;
                std::cout << "(iter,plane_id,cos_1,cos_2): " << iter << "," << plane_id << "," << cos_1 << "," << cos_2 << std::endl;

                is_parallel = true;
                std::get<0>(result) = is_parallel;
                return result;
            }
            factor = - refract_ratio * cos_1 - std::sqrt(cos_2);

            // get new line
            line.pt = pt_cross;
            line.unit_vector = line.unit_vector * refract_ratio + plane.norm_vector * factor;

            // get new plane
            plane.pt -=  plane.norm_vector * _pinplate_param.w_array[plane_id-1];

            // calculate new intersection
            is_parallel = myMATH::crossPoint(pt_cross, line, plane);
            if (is_parallel)
            {
                std::get<0>(result) = is_parallel;
                return result;
            }
        }
        pt_last = pt_cross;
    
        // backward projection: find intersection with the last plane
        plane_id = _pinplate_param.n_plate;

        // get line vector
        refract_ratio = _pinplate_param.refract_array[plane_id] / _pinplate_param.refract_array[plane_id+1];
        cos_1 = myMATH::dot(line.unit_vector, plane.norm_vector);
        cos_2 = 1 - refract_ratio*refract_ratio*(1 - cos_1*cos_1);
        if (cos_2 <= 0)
        {
            std::cout << "Camera::RefractPlate line " << __LINE__ << " : Error: total internal reflection" << std::endl;
            std::cout << "(iter,plane_id,cos_1,cos_2,delta): " << iter << "," << plane_id << "," << cos_1 << "," << cos_2 << "," << delta << std::endl;

            is_parallel = true;
            std::get<0>(result) = is_parallel;
            return result;
        }
        factor = - refract_ratio * cos_1 - std::sqrt(cos_2);

        line.unit_vector = line.unit_vector * refract_ratio + plane.norm_vector * factor;
        line.pt = _pinplate_param.t_vec_inv;

        // calculate intersection
        is_parallel = myMATH::crossPoint(pt_cross, line, plane);
        if (is_parallel)
        {
            std::get<0>(result) = is_parallel;
            return result;
        }

        // calculate diff vector
        delta = 0;
        for (int i = 0; i < 3; i ++)
        {
            diff_vec[i] = pt_cross[i] - pt_last[i];
            delta += diff_vec[i] * diff_vec[i];
        }
        
        // check convergence
        if (delta < _pinplate_param.proj_tol2)
        {
            std::get<1>(result) = pt_cross;
            std::get<2>(result) = delta;
            return result;
        }

        // update pt_init
        for (int i = 0; i < 3; i ++)
        {
            pt_init[i] += _pinplate_param.lr * diff_vec[i];
        }
    }

    std::get<0>(result) = is_parallel;
    std::get<1>(result) = pt_cross;
    std::get<2>(result) = delta;
    return result;
}

Pt2D Camera::distort (Pt2D const& pt_img_undist, PinholeParam const& param) const
{
    // opencv distortion model
    double x = pt_img_undist[0];
    double y = pt_img_undist[1];
    double xd, yd;

    // judge if distortion is needed
    if (param.is_distorted)
    {
        double r2 = x*x + y*y;
        double r4 = r2*r2;
        double r6 = r4*r2;

        double a1 = 2 * x*y;
        double a2 = r2 + 2 * x*x;
        double a3 = r2 + 2 * y*y;

        double cdist = 1 + param.dist_coeff[0]*r2 + param.dist_coeff[1]*r4;
        if (param.n_dist_coeff > 4)
        {
            cdist += param.dist_coeff[4]*r6;
        }         

        double icdist2 = 1.0;
        if (param.n_dist_coeff > 5)
        {
            icdist2 = 1.0 / (1 + param.dist_coeff[5]*r2 + param.dist_coeff[6]*r4 + param.dist_coeff[7]*r6);
        }
        // double icdist2 = 1.0 / (1 + param.dist_coeff[5]*r2 + param.dist_coeff[6]*r4 + param.dist_coeff[7]*r6);

        xd = x*cdist*icdist2 + param.dist_coeff[2]*a1 + param.dist_coeff[3]*a2;
        yd = y*cdist*icdist2 + param.dist_coeff[2]*a3 + param.dist_coeff[3]*a1;
        if (param.n_dist_coeff > 8)
        {
            // additional distortion by projecting onto a tilt plane (14)
            // not supported yet
            xd += param.dist_coeff[8] *r2 + param.dist_coeff[9] *r4;
            yd += param.dist_coeff[10]*r2 + param.dist_coeff[11]*r4;
        }
    }
    else
    {
        xd = x;
        yd = y;
    }
    
    Pt2D pt_img_dist(
        xd * param.cam_mtx(0,0) + param.cam_mtx(0,2),
        yd * param.cam_mtx(1,1) + param.cam_mtx(1,2)
    );

    return pt_img_dist;
}

// Polynomial model
Pt2D Camera::polyProject (Pt3D const& pt_world) const
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
Line3D Camera::lineOfSight (Pt2D const& pt_img_dist) const
{
    if (_type == PINHOLE)
    {
        return pinholeLine(undistort(pt_img_dist, _pinhole_param));
    }
    else if (_type == POLYNOMIAL)
    {
        return polyLineOfSight(pt_img_dist);
    }
    else if (_type == PINPLATE)
    {
        return pinplateLine(undistort(pt_img_dist, _pinplate_param));
    }
    else
    {
        std::cerr << "Camera::LineOfSight line " << __LINE__ << " : Error: unknown camera type: " << _type << std::endl;
        throw error_type;
    }
}

// Pinhole  model
Pt2D Camera::undistort (Pt2D const& pt_img_dist, PinholeParam const& param) const
{
    double fx = param.cam_mtx(0,0);
    double fy = param.cam_mtx(1,1);
    double cx = param.cam_mtx(0,2);
    double cy = param.cam_mtx(1,2);

    double u = pt_img_dist[0];
    double v = pt_img_dist[1];
    double x = (u - cx) / fx;
    double y = (v - cy) / fy;
    Pt2D img_undist(x,y);
    Pt2D img_dist(0,0);

    if (param.is_distorted)
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
            if (param.n_dist_coeff == 4)
            {
                icdist = 1 / (1 + (param.dist_coeff[1]*r2 + param.dist_coeff[0])*r2);
            }
            else if (param.n_dist_coeff == 5)
            {
                icdist = 1 / (1 + ((param.dist_coeff[4]*r2 + param.dist_coeff[1])*r2 + param.dist_coeff[0])*r2);
            }
            else
            {
                icdist = (1 + ((param.dist_coeff[7]*r2 + param.dist_coeff[6])*r2 + param.dist_coeff[5])*r2)/(1 + ((param.dist_coeff[4]*r2 + param.dist_coeff[1])*r2 + param.dist_coeff[0])*r2);
            }

            if (icdist < 0)
            {
                img_undist[0] = (u - cx) / fx;
                img_undist[1] = (v - cy) / fy;
                break;
            }

            // update x, y
            double deltaX = param.dist_coeff[2] * a1 + param.dist_coeff[3] * a2;
            double deltaY = param.dist_coeff[2] * a3 + param.dist_coeff[3] * a1;
            if (param.n_dist_coeff > 8)
            {
                deltaX += param.dist_coeff[8] * r2 + param.dist_coeff[9] * r2 * r2;
                deltaY += param.dist_coeff[10] * r2 + param.dist_coeff[11] * r2 * r2;
            }
            x = (x0 - deltaX) * icdist;
            y = (y0 - deltaY) * icdist;

            // calculate distorted xd, yd
            img_undist[0] = x;
            img_undist[1] = y;
            img_dist = distort(img_undist, param);

            // update error
            error = std::sqrt(std::pow(img_dist[0]-u, 2) + std::pow(img_dist[1]-v, 2));

            iter ++;
        }
    }

    return img_undist;
}

Line3D Camera::pinholeLine (Pt2D const& pt_img_undist) const
{
    Pt3D pt_world(pt_img_undist[0], pt_img_undist[1], 1.0);

    // calculate line of sight
    pt_world = _pinhole_param.r_mtx_inv * pt_world + _pinhole_param.t_vec_inv;

    Line3D line;
    line.pt = _pinhole_param.t_vec_inv;
    line.unit_vector = myMATH::createUnitVector(_pinhole_param.t_vec_inv, pt_world);

    return line;
}

Line3D Camera::pinplateLine (Pt2D const& pt_img_undist) const
{
    Pt3D pt_world(pt_img_undist[0], pt_img_undist[1], 1.0);

    // calculate line of sight
    pt_world = _pinplate_param.r_mtx_inv * pt_world + _pinplate_param.t_vec_inv;
    
    Line3D line;
    line.pt = _pinplate_param.t_vec_inv;
    line.unit_vector = myMATH::createUnitVector(_pinplate_param.t_vec_inv, pt_world);

    // backward projection
    Plane3D plane;
    plane.pt = _pinplate_param.pt3d_closest;
    plane.norm_vector = _pinplate_param.plane.norm_vector;

    Pt3D pt_cross;
    bool is_parallel = false;
    double refract_ratio;
    double cos_1;
    double cos_2;
    double factor;

    for (int plane_id = _pinplate_param.n_plate; plane_id >= 0; plane_id --)
    {
        is_parallel = myMATH::crossPoint(pt_cross, line, plane);
        if (is_parallel)
        {
            std::cerr << "Camera::PinplateLine line " << __LINE__ << " : Error: line is parallel to the plane" << std::endl;
            throw error_type;
        }

        // if (plane_id == _pinplate_param.n_plate)
        // {
        //     line.unit_vector = myMATH::createUnitVector(pt_cross, _pinplate_param.t_vec_inv);
        // }

        // get new line
        refract_ratio = _pinplate_param.refract_array[plane_id+1] / _pinplate_param.refract_array[plane_id];
        cos_1 = myMATH::dot(line.unit_vector, plane.norm_vector);
        cos_2 = 1 - refract_ratio*refract_ratio*(1 - cos_1*cos_1);
        if (cos_2 <= 0)
        {
            std::cerr << "Camera::PinplateLine line " << __LINE__ << " : Error: total internal reflection" << std::endl;
            throw error_type;
        }
        factor = - refract_ratio * cos_1 + std::sqrt(cos_2);

        line.unit_vector = line.unit_vector * refract_ratio + plane.norm_vector * factor;
        line.pt = pt_cross;

        // update plane location
        if (plane_id > 0)
        {
            plane.pt += plane.norm_vector * _pinplate_param.w_array[plane_id-1];
        } 
    }

    return line;
}

// Polynomial model
Pt3D Camera::polyImgToWorld (Pt2D const& pt_img_dist, double plane_world) const
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
        // pt_img_temp = polyProject(pt_world);
        // du = pt_img_dist[0] - pt_img_temp[0];
        // dv = pt_img_dist[1] - pt_img_temp[1];

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

Line3D Camera::polyLineOfSight (Pt2D const& pt_img_dist) const
{
    Pt3D pt_world_1 = polyImgToWorld(pt_img_dist, _poly_param.plane[0]);
    Pt3D pt_world_2 = polyImgToWorld(pt_img_dist, _poly_param.plane[1]);
    Pt3D unit_vec = myMATH::createUnitVector(pt_world_1, pt_world_2);

    Line3D line = {pt_world_1, unit_vec};

    return line;
}