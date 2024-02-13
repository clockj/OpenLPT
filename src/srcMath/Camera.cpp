#include "Camera.h"

Camera::Camera () {};

Camera::Camera (const Camera& c)
    : _n_off_h(c._n_off_h), _n_off_w(c._n_off_w), _n_pix_w(c._n_pix_w), 
      _n_pix_h(c._n_pix_h), _w_pix(c._w_pix), _h_pix(c._h_pix), 
      _f_eff(c._f_eff), _kr(c._kr), _kx(c._kx), _r_mtx(c._r_mtx), 
      _t_vec(c._t_vec), _r_mtx_inv(c._r_mtx_inv), _t_vec_inv(c._t_vec_inv) {}

Camera::Camera(std::istream& is)
    : _r_mtx(3), _r_mtx_inv(3), _t_vec(3,1), _t_vec_inv(3,1)
{
    LoadParameters (is);
}

Camera::Camera(std::string file_name)
    : _r_mtx(3), _r_mtx_inv(3), _t_vec(3,1), _t_vec_inv(3,1)
{
    LoadParametersFromFile(file_name);
}

void Camera::LoadParameters (std::istream& is)
{
    is >> _n_off_h;
    is >> _n_off_w;
    is >> _n_pix_w;
    is >> _n_pix_h;
    is >> _w_pix;
    is >> _h_pix;
    is >> _f_eff;
    is >> _kr;
    is >> _kx;

    _kr = std::fabs(_kr);
    _kx = std::fabs(_kx);

    double val = 0.0;
    for (int i = 0; i < 3; i ++)
    {
        for (int j = 0; j < 3; j ++)
        {
            is >> val;
            _r_mtx(i,j) = val;
        }
    }
    for (int i = 0; i < 3; i ++) {
        is >> val;
        _t_vec(i,0) = val;
    }
    for (int i = 0; i < 3; i ++)
    {
        for (int j = 0; j < 3; j ++)
        {
            is >> val;
            _r_mtx_inv(i,j) = val;
        }
    }
    for (int i = 0; i < 3; i ++) {
        is >> val;
        _t_vec_inv(i,0) = val;
    }
}

void Camera::LoadParametersFromFile (std::string file_name)
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
            file_content << line << '\t';
    	}
    }
    infile.close();

    LoadParameters(file_content);
}

//
// Projection
//

Matrix<double> Camera::ImgPixelToMM (Matrix<double> const& pt_img_pix)
{
    // Position centered (p);
    //  shift origin to the center of the image
    //  centered -= Position(Npixw/2, Npixh/2, 0);

    // account for left-handed coordinate system 
    Matrix<double> pt_img_mm (pt_img_pix); 
    // Tsai's model: xf,yf start from 1
    // But here, xf,yf start from 0
    pt_img_mm(0,0) += 1;
    pt_img_mm(1,0) += 1;
    pt_img_mm(0,0) = _kx*_w_pix * (  pt_img_mm(0,0) - _n_pix_w/2 - _n_off_w);
    pt_img_mm(1,0) = _h_pix     * (- pt_img_mm(1,0) + _n_pix_h/2 - _n_off_h);
    pt_img_mm(2,0) = 0;

    if (_kr > SMALLNUMBER) 
    {
        // Tsai's model 
        // (the original code is wrong) 
        double x = pt_img_mm(0,0);
        double y = pt_img_mm(1,0);
        double r_square = x*x + y*y;
    
        x = x * (1 + _kr * r_square); // Xu=Xd+Dx, Dx=Xd(1+k_r*r^2), r=sqrt(Xd^2+Yd^2)
        y = y * (1 + _kr * r_square); // Yu=Yd+Dy, Dy=Yd(1+k_r*r^2), r=sqrt(Xd^2+Yd^2)
        
        pt_img_mm(0,0) = x;
        pt_img_mm(1,0) = y;
        return pt_img_mm;
    }
    else
    {
        return pt_img_mm;
    }
}

Matrix<double> Camera::ImgMMToPixel (Matrix<double> const& pt_img_mm)
{
    // Tsai's model 
    // Xu=Xd+Dx, Dx=Xd(1+k_r*r^2), r=sqrt(Xd^2+Yd^2)
    // Yu=Yd+Dy, Dy=Yd(1+k_r*r^2), r=sqrt(Xd^2+Yd^2) 
    // r_u_2 = xu^2 + yu^2
    // compute Xd
    // Xd_1_num = - (2/3)^(1/3) * Xu^2 
    // Xd_1_den = ( 9 * k1^2 * Xu^3 * (r_u_2)^2 + sqrt(3) * sqrt( k1^3*Xu^6(r_u_2)^3 * (4 + 27 * k1 * r_u_2) ) )^(1/3)
    // Xd_2_num = Xd_1_den
    // Xd_2_den = 2^(1/3) * 3^(2/3) * k1 * r_u_2
    // Xd = Xd_1_num / Xd_1_den + Xd_2_num / Xd_2_den
    // Compute Yd 
    // Yd_1_num = - (2/3)^(1/3) * Xu * Yu 
    // Yd_1_den = Xd_1_den
    // Yd_2_num = Yu * Xd_1_den
    // Yd_2_den = Xu * Xd_2_den
    // Yd = Yd_1_num / Yd_1_den + Yd_2_num / Yd_2_den 
    // Note: 
    // we need to make extra cases for Xu = 0 or Yu = 0

    double x_u = pt_img_mm[0];
    double y_u = pt_img_mm[1];

    double x_d = x_u;
    double y_d = y_u;

    // remove potential cylindrical distortion
    if (_kr > SMALLNUMBER)
    {
        // for solve x_d
        double x_d_1_num = 0.0;
        double x_d_1_den = 0.0;
        double x_d_2_num = 0.0;
        double x_d_2_den = 0.0;
        // for solve y_d
        double y_d_1_num = 0.0;
        double y_d_1_den = 0.0;
        double y_d_2_num = 0.0;
        double y_d_2_den = 0.0;

        if (std::fabs(x_u) < SMALLNUMBER || std::fabs(y_u) < SMALLNUMBER)
        {
            // x_u is very small 
            // or both x_u and y_u are very small 
            if (std::fabs(x_u) < SMALLNUMBER)
            {
                x_d = 0.0;
                if (std::fabs(y_u) < SMALLNUMBER)
                {
                    // both x_u, y_u are small 
                    y_d = 0.0;
                }
                else 
                {
                    // x_u is small
                    // y_u is large enough
                    y_d_1_num = - std::pow(2./3., 1./3.);
                    y_d_1_den = std::pow( 9 * _kr*_kr * y_u + std::sqrt(3) * std::sqrt( 4 * std::pow(_kr, 3) + 27 * std::pow(_kr, 4) * y_u*y_u ), 1./3.);

                    y_d_2_num = y_d_1_den;
                    y_d_2_den = std::pow(2, 1./3.) * std::pow(3, 2./3.) * _kr;

                    y_d = y_d_1_num / y_d_1_den + y_d_2_num / y_d_2_den;
                }
            }
            else 
            {
                // x_u is large enough 
                // y_u is small 
                y_d = 0.0;

                x_d_1_num = - std::pow(2./3., 1./3.);
                x_d_1_den = std::pow( 9 * _kr*_kr * x_u + std::sqrt(3) * std::sqrt( 4 * std::pow(_kr, 3) + 27 * std::pow(_kr, 4) * x_u*x_u ), 1./3.);

                x_d_2_num = x_d_1_den;
                x_d_2_den = std::pow(2, 1./3.) * std::pow(3, 2./3.) * _kr;

                x_d = x_d_1_num / x_d_1_den + x_d_2_num / x_d_2_den;
            }
        }
        else
        {
            double r_u_2 = x_u * x_u + y_u * y_u;

            x_d_1_num = std::pow(2./3., 1./3.) * x_u*x_u;
            x_d_1_den = std::pow( 
                            9 * _kr*_kr * std::pow(x_u,3) * r_u_2*r_u_2 
                            + std::sqrt(3)
                            * std::sqrt( 
                                std::pow(_kr,3)*std::pow(x_u,6)*std::pow(r_u_2,3) 
                                * (4 + 27 * _kr * r_u_2) 
                              ), 
                            1./3. );
            x_d_2_num = x_d_1_den;
            x_d_2_den = std::pow(2, 1./3.) * std::pow(3, 2./3.) * _kr * r_u_2;
            x_d = - x_d_1_num / x_d_1_den + x_d_2_num / x_d_2_den;

            y_d_1_num = std::pow(2./3., 1./3.) * x_u * y_u;
            y_d_1_den = x_d_1_den;
            y_d_2_num = y_u * x_d_2_num;
            y_d_2_den = x_u * x_d_2_den;
            y_d = - y_d_1_num / y_d_1_den + y_d_2_num / y_d_2_den;
        }
    }
    
    
    // remove Sx: uncertainty along x direction
    // scale into pixel units
    double xf = x_d/_kx / _w_pix;
    double yf = y_d     / _h_pix;

    // shift origin
    // origin is on the top left corner
    // ----> x direction (rightwards)
    // | 
    // |
    // \/ y direction (downwards)
    // Tsai's model: xf,yf start from 1
    xf = xf + _n_pix_w/2.0 + _n_off_w;
    yf = -1.0 * (yf - _n_pix_h/2.0) - _n_off_h;
    // But here, xf,yf start from 0
    xf -= 1;
    yf -= 1;
    return Matrix<double> ( 
        3,1,
        xf, yf, 0
    );
}

Matrix<double> Camera::ImgMMToWorld (Matrix<double> const& pt_img_mm)
{
    //Xu = f * x / z, Yu = f * y / z; (x, y, z): pt_cam (scaled), (Xu, Yu, 0): pt_img_mm (unscaled)
    Matrix<double> pt_cam(pt_img_mm);
    pt_cam *= _t_vec(2,0) / _f_eff;
    pt_cam(2,0) = _t_vec(2,0);

    // R^-1 * (pt_cam - T): pt_cam in world coordinate
    pt_cam -= _t_vec; 
    return _r_mtx_inv * pt_cam;
}

Matrix<double> Camera::WorldToImgMM (Matrix<double> const& pt_world)
{
    // TODO: add comments 
    Matrix<double> pt_img_mm = _r_mtx * pt_world + _t_vec;
    pt_img_mm *= _f_eff / pt_img_mm(2,0);
    pt_img_mm(2,0) = 0.0;
    return pt_img_mm;
}


