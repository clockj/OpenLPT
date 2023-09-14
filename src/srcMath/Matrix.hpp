#ifndef MATRIX_HPP
#define MATRIX_HPP

#include "Matrix.h"


//
// Constructor
//
template<class T> 
Matrix<T>::Matrix(int dim_x, int dim_y)
    : Matrix(dim_x, dim_y, 0)
{}

template<class T> 
Matrix<T>::Matrix(int dim_x, int dim_y, T val)
    : _dim_x(dim_x), _dim_y(dim_y)
{
    CreateMtx(val);
}

template<class T> 
Matrix<T>::Matrix (std::vector<std::vector<T>> const& mtx)
    : _dim_x(mtx.size()), _dim_y(mtx[0].size())
{
    _n = _dim_x * _dim_y;
    _is_space = 1;
    _mtx = new T [_n];

    for (int i = 0; i < _dim_x; i ++)
    {
        for (int j = 0; j < _dim_y; j ++)
        {
            _mtx[MapID(i,j)] = mtx[i][j];
        }
    }
}

template<class T> 
Matrix<T>::Matrix (const Matrix<T>& mtx) 
    : _dim_x(mtx._dim_x), _dim_y(mtx._dim_y), _n(mtx._n), _is_space(1)
{
    _mtx = new T [_n];
    for (int i = 0; i < _n; i ++)
    {
        _mtx[i] = mtx._mtx[i];
    }
}

// template<class T> 
// Matrix<T>::Matrix (int dim_x, int dim_y, const T* mtx) 
//     : _dim_x(dim_x), _dim_y(dim_y), _n(dim_x*dim_y)
// {
//     _mtx = new T [_n];
//     for (int i = 0; i < _n; i ++)
//     {
//         _mtx[i] = mtx[i];
//     }
// }

template<class T>
Matrix<T>::Matrix (const Matrix<T>& vec_1, const Matrix<T>& vec_2) // vec_1, vec_2 have same dim
    : _dim_x(vec_1._dim_x), _dim_y(vec_1._dim_y), _n(vec_1._n), _is_space(1)
{
    if (vec_2._dim_x != _dim_x || vec_2._dim_y != _dim_y)
    {
        std::cerr << "Matrix::Matrix: the sizes of vec_1 and vec_2 are different" 
                  << std::endl;
        throw error_size;
    }

    _mtx = new T [_n];

    T mag = 0;
    T value = 0;
    for (int i = 0; i < _n; i ++)
    {
        value = vec_1._mtx[i] - vec_2._mtx[i];
        _mtx[i] = value;
        mag += (value * value);
    }

    if (mag < MAGSMALLNUMBER)
    {
        std::cerr << "Matrix::Matrix: vec_1 and vec_2 are too close"
                  << std::endl;
        throw error_div0;
    }

    mag = T (std::sqrt(mag));
    for (int i = 0; i < _n; i ++)
    {
        _mtx[i] /= mag;
    }
}

template<class T> 
Matrix<T>::Matrix (T m11, T m12,
                   T m21, T m22)
    : _dim_x(2), _dim_y(2), _n(4), _is_space(1)
{
    _mtx = new T [4];
    _mtx[0] = m11; _mtx[1] = m12;
    _mtx[2] = m21; _mtx[3] = m22;
}

template<class T> 
Matrix<T>::Matrix (T m11, T m12, T m13,
                   T m21, T m22, T m23,
                   T m31, T m32, T m33)
    : _dim_x(3), _dim_y(3), _n(9), _is_space(1)
{    
    _mtx = new T [9];

    _mtx[0] = m11; _mtx[1] = m12; _mtx[2] = m13;
    _mtx[3] = m21; _mtx[4] = m22; _mtx[5] = m23;
    _mtx[6] = m31; _mtx[7] = m32; _mtx[8] = m33;
}

template<class T>
Matrix<T>::Matrix (int dim_x, int dim_y, T v1, T v2, T v3)
    : _dim_x(3), _dim_y(1), _n(3), _is_space(1)
{
    if (dim_x != 3 || dim_y != 1)
    {
        std::cerr << "The dimension does not satisfied for "
                  << "creating a 3*1 vector!"
                  << std::endl;
        throw error_size;
    }
    _mtx = new T [3];

    _mtx[0] = v1;
    _mtx[1] = v2;
    _mtx[2] = v3;
}

template<class T> 
Matrix<T>::Matrix (int dim)
    : Matrix(dim, dim, 0)
{
    for (int i = 0; i < _dim_x; i ++)
    {
        _mtx[MapID(i,i)] = 1;
    }
}

template<class T>
Matrix<T>::Matrix (std::string file_name)
{
    // std::cout << "Start loading!" << std::endl;

    std::string line;
    T value;

    // 1st loop: get _dim_x, _dim_y
    int dim_x = 0;
    int dim_y = 0;
    
    std::ifstream infile;
    std::istringstream istream;
    infile.open(file_name);
    while (std::getline(infile, line))
    {       
        dim_x ++;
        
        if (dim_x == 1)
        {
            istream.str(line);
            while(istream >> value)
            {
                dim_y ++;
                if (istream.peek() == ',')
                {
                    istream.ignore();
                }
            }
        }
    }
    infile.close();
    

    if (dim_x!=0 && dim_y!=0) 
    {
        CreateMtx(dim_x, dim_y, 0);
        
        // 2nd loop: load data
        std::ifstream infile_data;
        infile_data.open(file_name);
        
        for (int i = 0; i < _dim_x; i ++)
        {
            std::getline(infile_data, line);
            std::istringstream istream_data;
            istream_data.str(line);

            int id_start = i * _dim_y;
            int j = 0;
            while (istream_data >> value)
            {
                _mtx[id_start+j] = value;
                j ++;
                if (istream_data.peek() == ',')
                {
                    istream_data.ignore();
                }
            } // NOTE: make sure each row has the same elem
        }

        infile_data.close();
    }
    else 
    {
        CreateMtx(1,1,0);
    }

    // std::cout << "Finish loading!" << std::endl;
}

// Deconstructor 
template<class T>
Matrix<T>::~Matrix ()
{
    ClearMtx();
}


// Create/clear space 
template<class T>
void Matrix<T>::ClearMtx()
{
    if (_is_space)
    {
        _dim_x = 0;
        _dim_y = 0;
        _n = 0;
        _is_space = 0;

        delete[] _mtx;
    }
}

template<class T>
void Matrix<T>::CreateMtx(int dim_x, int dim_y, T val)
{
    if (!_is_space)
    {
        _is_space = 1;

        _dim_x = dim_x;
        _dim_y = dim_y;
        _n = _dim_x * _dim_y;

        _mtx = new T [_n];

        for (int i = 0; i < _n; i ++)
        {
            _mtx[i] = val;
        }
    }
    else
    {
        std::cerr << "The space for _mtx is not cleared: "
                  << "Cannot claim more space!" << std::endl;
        throw error_space;
    }
}

template<class T>
void Matrix<T>::CreateMtx(T val)
{
    CreateMtx(_dim_x, _dim_y, val);
}


// Get id
template<class T>
inline int Matrix<T>::MapID(int id_x, int id_y) const
{
    return (id_x * _dim_y + id_y);
}


// Type transform
template<class T> 
Matrix<double> Matrix<T>::TypeToDouble ()
{
    Matrix<double> res(_dim_x,_dim_y);
    for (int i = 0; i < _n; i ++)
    {
        res._mtx[i] = double(_mtx[i]);
    }
    return res;
}

template<class T> 
Matrix<int> Matrix<T>::TypeToInt ()
{
    Matrix<int> res(_dim_x, _dim_y);
    for (int i = 0; i < _n; i ++)
    {
        res._mtx[i] = int(_mtx[i]);
    }
    return res;
}



//
// Get/Assign value
//
template<class T> 
T Matrix<T>::operator() (int i, int j) const
{
    return _mtx[MapID(i,j)];
}

template<class T> 
T& Matrix<T>::operator() (int i, int j)
{
    return _mtx[MapID(i,j)];
}

// Just for vector
template<class T>
T Matrix<T>::operator[] (int i) const 
{
    return _mtx[i];
}

// Just for vector
template<class T>
T& Matrix<T>::operator[] (int i)
{
    return _mtx[i];
}

template<class T> 
std::vector<T> Matrix<T>::GetRow(int i)
{
    std::vector<T> row(_dim_y, 0);
    int id_start = i * _dim_y;
    for (int iter = 0; iter < _dim_y; iter ++)
    {
        row[iter] = _mtx[id_start + iter];
    }
    return row;
}

template<class T> 
std::vector<T> Matrix<T>::GetCol(int j)
{
    std::vector<T> col(_dim_x, 0);
    for (int iter = 0; iter < _dim_x; iter ++)
    {
        col[iter] = _mtx[MapID(iter,j)];
    }
    return col;
}


//
// Get matrix info
//
template<class T> 
inline int Matrix<T>::GetDimX() const
{
    return _dim_x;
}

template<class T> 
inline int Matrix<T>::GetDimY() const
{
    return _dim_y;
}

template<class T> 
void Matrix<T>::Print()
{
    std::cout << std::endl;
    std::cout << "Matrix = " << std::endl;
    std::cout << "(dim_x, dim_y) = " 
              << "(" << _dim_x << "," << _dim_y << ")" 
              << std::endl;
    int id_start;
    for (int iter_x = 0; iter_x < _dim_x; iter_x ++)
    {
        id_start = iter_x * _dim_y;

        std::cout << "| ";
        for (int iter_y = 0; iter_y < _dim_y; iter_y ++)
        {
            std::cout << std::scientific 
                      << std::setprecision(2)  
                      << std::setfill(' ') 
                      << std::left 
                      << _mtx[id_start + iter_y] << " ";
        }
        std::cout << "|" << std::endl;
    }
    std::cout << std::endl;
}

template<class T> 
double Matrix<T>::Mag()
{
    return FrobeniusNorm();
}

template<class T> 
double Matrix<T>::FrobeniusNorm()
{
    double res = 0.0;
    for (int i = 0; i < _n; i ++)
    {
        res += (_mtx[i] * _mtx[i]);
    }
    return std::sqrt(res);
}

template<class T> 
double Matrix<T>::Dot(Matrix<T> const& mtx)
{
    if (_dim_x != mtx._dim_x || _dim_y != mtx._dim_y)
    {
        std::cerr << "The dimension does not match: "
                  << "cannot calculate dot!" 
                  << std::endl;
        throw error_size;
    }

    double res = 0;
    for (int i = 0; i < _n; i ++)
    {
        res += (_mtx[i] * mtx._mtx[i]);
    }

    return res;
}

template<class T>
double Matrix<T>::DistSqr(Matrix<T> const& mtx)
{
    if (_dim_x!=mtx._dim_x || _dim_y!=mtx._dim_y)
    {
        std::cerr << "Matrix::Dist: size not matched!" << std::endl;
        throw error_size;
    }

    double dist = 0;
    for (int i = 0; i < _n; i ++)
    {
        dist += std::pow((double)_mtx[i] - mtx._mtx[i], 2);
    }

    return dist;
}

template<class T>
double Matrix<T>::DistSqr(Matrix<T> const& pt, Matrix<T> const& unit)
{
    if (_dim_x!=3||_dim_y!=1||pt._dim_x!=3||pt._dim_y!=1||unit._dim_x!=3||unit._dim_y!=1)
    {
        std::cerr << "Matrix::Dist: size not matched!" << std::endl;
        throw error_size;
    }

    double x1 = pt._mtx[0] - _mtx[0];
    double x2 = pt._mtx[1] - _mtx[1];
    double x3 = pt._mtx[2] - _mtx[2];

    // calculate dot 
    double res = x1*unit._mtx[0] + x2*unit._mtx[1] + x3*unit._mtx[2];

    double dist = x1*x1 + x2*x2 + x3*x3 - res*res;

    if (dist < 0 && dist >= - MAGSMALLNUMBER)
    {
        dist = 0.0;
    }
    else if (dist < - MAGSMALLNUMBER)
    {
        std::cout << " error: "
            << "the distance in  is: "
            << "dist = "
            << dist
            << std::endl;
        throw error_range;
    }

    return dist;
}

template<class T>
double Matrix<T>::Dist(Matrix<T> const& mtx)
{
    return std::sqrt(DistSqr(mtx));
}

template<class T>
double Matrix<T>::Dist(Matrix<T> const& pt, Matrix<T> const& unit)
{
    return std::sqrt(DistSqr(pt, unit));
}


//
// Matrix output 
//
template<class T>
void Matrix<T>::WriteMatrix (std::string file_name)
{
    std::cout << "\nStart writing!" << std::endl;

    std::ofstream outfile(file_name, std::ios::out);

    outfile << "Matrix: " << "(row,col)=" 
            << "(" << _dim_x << "," 
            << _dim_y << ")" << "\n";

    T value;
    for (int i = 0; i < _dim_x; i ++)
    {
        int id_start = i * _dim_y;
        for (int j = 0; j < _dim_y-1; j ++)
        {
            value = _mtx[id_start + j];
            outfile << value << ",";
        }
        outfile << _mtx[id_start + _dim_y-1] << "\n";
    }

    outfile.close();
    std::cout << "Finish writing!" << std::endl;
}


//
// Matrix calculation
//

template<class T> 
Matrix<T>& Matrix<T>::operator= (Matrix<T> const& mtx)
{
    if (_dim_x != mtx._dim_x || _dim_y != mtx._dim_y)
    {
        ClearMtx();

        _is_space = 1;
        _dim_x = mtx._dim_x;
        _dim_y = mtx._dim_y;
        _n = mtx._n;

        _mtx = new T [_n];
    }

    for (int i = 0; i < _n; i ++)
    {
        _mtx[i]= mtx._mtx[i];
    }

    return *this;
}

template<class T> 
Matrix<T> Matrix<T>::operator+ (Matrix<T> const& mtx)
{
    if ( (_dim_x != mtx._dim_x) || (_dim_y != mtx._dim_y) )
    {
        std::cerr << "The size of matrices do not match!" << std::endl;
        throw error_size;
    }

    Matrix<T> res(mtx);
    for (int i = 0; i < _n; i ++)
    {
        res._mtx[i] += _mtx[i];
    }

    return res;
}

template<class T> 
Matrix<T>& Matrix<T>::operator+= (Matrix<T> const& mtx)
{
    if ((_dim_x != mtx._dim_x) || (_dim_y != mtx._dim_y))
    {
        std::cerr << "The size of matrices do not match!" << std::endl;
        throw error_size;
    }
    for (int i = 0; i < _n; i ++)
    {
        _mtx[i] += mtx._mtx[i];
    }
    return *this;
}

template<class T> 
Matrix<T> Matrix<T>::operator- (Matrix<T> const& mtx)
{
    if ( (_dim_x != mtx._dim_x) || (_dim_y != mtx._dim_y) )
    {
        std::cerr << "The size of matrices do not match!" << std::endl;
        throw error_size;
    }

    Matrix<T> res(_dim_x, _dim_y, 0);

    for (int i = 0; i < _n; i ++)
    {
        res._mtx[i] = _mtx[i] - mtx._mtx[i];
    }

    return res;
}

template<class T> 
Matrix<T>& Matrix<T>::operator-= (Matrix<T> const& mtx)
{
    if ((_dim_x != mtx._dim_x) || (_dim_y != mtx._dim_y))
    {
        std::cerr << "The size of matrices do not match!" << std::endl;
        throw error_size;
    }
    for (int i = 0; i < _n; i ++)
    {
        _mtx[i] -= mtx._mtx[i];
    }
    return *this;
}

template<class T> 
Matrix<T> Matrix<T>::operator* (Matrix<T> const& mtx)
{
    if ( _dim_y != mtx._dim_x )
    {
        std::cerr << "The size of matrices do not match!" << std::endl;
        throw error_size;
    }

    Matrix<T> res(_dim_x, mtx._dim_y, 0);
    for (int iter_x = 0; iter_x < res._dim_x; iter_x ++)
    {
        for (int iter_y = 0; iter_y < res._dim_y; iter_y ++)
        {
            for (int iter = 0; iter < _dim_y; iter ++)
            {
                res(iter_x,iter_y) += (_mtx[MapID(iter_x,iter)] * mtx(iter,iter_y));
            }
        }
    }

    return res;
}

template<class T> 
Matrix<T>& Matrix<T>::operator*= (Matrix<T> const& mtx)
{
    if ( _dim_y != mtx._dim_x )
    {
        std::cerr << "The size of matrices do not match!" << std::endl;
        throw error_size;
    }

    Matrix<T> res(_dim_x, mtx._dim_y, 0);
    for (int iter_x = 0; iter_x < res._dim_x; iter_x ++)
    {
        for (int iter_y = 0; iter_y < res._dim_y; iter_y ++)
        {
            for (int iter = 0; iter < _dim_y; iter ++)
            {
                res(iter_x,iter_y) += (_mtx[MapID(iter_x,iter)] * mtx(iter,iter_y));
            }
        }
    }
    
    if (_dim_y != res._dim_y)
    {
        // _dim_x = res._dim_x;
        _dim_y = res._dim_y;
        _n = res._n;

        delete[] _mtx;
        _mtx = new T [_n];
    }

    for (int i = 0; i < _n; i ++)
    {
        _mtx[i] = res._mtx[i];
    }

    return *this;
}


template<class T> 
Matrix<T> Matrix<T>::operator* (T ratio)
{
    Matrix<T> res(*this);
    for (int i = 0; i < _n; i ++)
    {
        res._mtx[i] *= ratio;
    }
    return res;
}

template<class T> 
Matrix<T>& Matrix<T>::operator*= (T ratio)
{
    for (int i = 0; i < _n; i ++)
    {
        _mtx[i] *= ratio;   
    }
    return *this;
}

template<class T> 
Matrix<T> Matrix<T>::operator/ (T ratio)
{
    Matrix<T> res(*this);
    for (int i = 0; i < _n; i ++)
    {
        res._mtx[i] /= ratio;
    }
    return res;
}

template<class T> 
Matrix<T>& Matrix<T>::operator/= (T ratio)
{
    for (int i = 0; i < _n; i ++)
    {
        _mtx[i] /= ratio;
    }
    return *this;
}


template<class T> 
void Matrix<T>::PiecewiseMultiply (Matrix<T> const& mtx)
{
    // TODO: justify dimensions
    for (int i = 0; i < _n; i ++)
    {
        _mtx[i] *= mtx._mtx[i];
    }
}



//
// Matrix manipulation
//
template<class T> 
Matrix<T> Matrix<T>::Transpose ()
{
    Matrix<T> res(_dim_y, _dim_x, 0);
    for (int iter_x = 0; iter_x < _dim_x; iter_x ++)
    {
        for (int iter_y = 0; iter_y < _dim_y; iter_y ++)
        {
            res(iter_y,iter_x) = _mtx[MapID(iter_x,iter_y)];
        }
    }
    return res;
}

template<class T> 
Matrix<T> Matrix<T>::Diagonal ()
{
    if (_dim_x != _dim_y)
    {
        std::cerr << "It is not a square matrix;" 
                  << "there is no diagonal matrix!" 
                  << std::endl;
        throw error_size;
    }

    Matrix<T> res(_dim_x, _dim_x, 0);
    for (int iter_x = 0; iter_x < _dim_x; iter_x ++)
    {
        res(iter_x,iter_x) = _mtx[MapID(iter_x,iter_x)];
    }

    return res;
}

template<class T> 
double Matrix<T>::Trace ()
{
    if (_dim_x != _dim_y)
    {
        std::cerr << "It is not a square matrix;" 
                  << "there is no trace!" 
                  << std::endl;
        throw error_size;
    }

    double res = 0.0;
    for (int iter_x = 0; iter_x < _dim_x; iter_x ++)
    {
        res += _mtx[MapID(iter_x,iter_x)];
    }

    return res;
}

template<class T> 
Matrix<double> Matrix<T>::Inverse ()
{
    return GaussInverse();
}

template<class T> 
Matrix<double> Matrix<T>::GaussInverse ()
{
    if (_dim_x != _dim_y)
    {
        std::cerr << "The matrix is not a square matrix!"
                  << "There is no inverse."
                  << std::endl;
        throw error_size;
    }

    // Gaussian Elimination
    std::vector<int> sortedRowIndex(_dim_x);
    Matrix<double> u_mtx = this->TypeToDouble(); // deep copy
    Matrix<double> b_mtx(_dim_x); // generate identity matrix
    Matrix<double> res(_dim_x);

    GaussForward(u_mtx, sortedRowIndex);
    res = GaussBackward(b_mtx, u_mtx, sortedRowIndex);

    return res;
}


template<class T> 
Matrix<double> Matrix<T>::DetInverse ()
{
    if (_dim_x != _dim_y)
    {
        std::cerr << "The matrix is not a square matrix!"
                  << "There is no inverse."
                  << std::endl;
        throw error_size;
    }

    if (_dim_x != 3)
    {
        std::cerr << "DetInverse only supports for 3x3 matrix!"
                  << std::endl;
        throw error_size;
    }
    
    Matrix<double> res(3,3,0);
    //        [0]    [1]    [2]
    // [0] res[0] res[1] res[2]
    // [1] res[3] res[4] res[5]
    // [2] res[6] res[7] res[8]

    double det =   _mtx[0] * ( _mtx[4]*_mtx[8] - _mtx[5]*_mtx[7] )
                 - _mtx[1] * ( _mtx[3]*_mtx[8] - _mtx[5]*_mtx[6] ) 
                 + _mtx[2] * ( _mtx[3]*_mtx[7] - _mtx[4]*_mtx[6] );
    res._mtx[0] =   (_mtx[4]*_mtx[8] - _mtx[5]*_mtx[7]) / det;
    res._mtx[1] = - (_mtx[3]*_mtx[8] - _mtx[5]*_mtx[6]) / det;
    res._mtx[2] =   (_mtx[3]*_mtx[7] - _mtx[4]*_mtx[6]) / det;

    res._mtx[3] = - (_mtx[1]*_mtx[8] - _mtx[2]*_mtx[7]) / det;
    res._mtx[4] =   (_mtx[0]*_mtx[8] - _mtx[2]*_mtx[6]) / det;
    res._mtx[5] = - (_mtx[0]*_mtx[7] - _mtx[1]*_mtx[6]) / det;

    res._mtx[6] =   (_mtx[1]*_mtx[5] - _mtx[2]*_mtx[4]) / det;
    res._mtx[7] = - (_mtx[0]*_mtx[5] - _mtx[2]*_mtx[3]) / det;
    res._mtx[8] =   (_mtx[0]*_mtx[4] - _mtx[1]*_mtx[3]) / det;

    return res;
}

template<class T> 
void Matrix<T>::GaussForward(Matrix<double>& u_mtx, std::vector<int>& sortedRowIndex)
{
    // Solve Amat * Xmat = Bmat,
    //  where Amat (n*n), Xmat (n*m), Bmat (n*m).
    // In GaussForward, Amat = Lmat * u_mtx, 
    //  we are trying to solve for u_mtx (upper triangle matrix).
    // The actions (Lmat^-1) are stored in 
    //  the lower triangle part of u_mtx (multipliers),
    //  and the sortedRowIndex.

    double maxPivot; // "biggest pivot" of each col
    int t;           // row index of "biggest pivot"
    double pivot;    // "relative" pivot value
    double m;        // multiplier during elimination

    std::vector<T> row(u_mtx._dim_y, 0.0);
    std::vector<T> maxOfrow(u_mtx._dim_x, 0.0);
    // std::vector<int> sortedRowIndex(_dim_x, 0);

    for (int iter_x = 0; iter_x < u_mtx._dim_x; iter_x ++)
    {
        sortedRowIndex[iter_x] = iter_x;
        row = GetRow(iter_x);
        maxOfrow[iter_x] = myMATH::Max<T>(row);
    }

    for (int iter_y = 0; iter_y < u_mtx._dim_y-1; iter_y ++)
    {
        maxPivot = 0.0; // initialize "biggest pivot" 
        t = 0;
        for (int iter_x = iter_y; iter_x < u_mtx._dim_x; iter_x ++)
        {
            pivot = std::fabs(u_mtx(sortedRowIndex[iter_x], iter_y) 
                              / maxOfrow[sortedRowIndex[iter_x]]);
            if (pivot > maxPivot)
            {
                maxPivot = pivot;
                t = iter_x;
            }
        }
        myMATH::Swap<int>(sortedRowIndex, t, iter_y); // move the row with "biggest pivot" 
                                                 // to the first row of the sub-matrix

        for (int iter_x = iter_y+1; iter_x < u_mtx._dim_x; iter_x ++)
        {
            m = u_mtx(sortedRowIndex[iter_x], iter_y) 
                / u_mtx(sortedRowIndex[iter_y], iter_y);
            u_mtx(sortedRowIndex[iter_x], iter_y) = m;  // it should be 0. 
                                                    // here is only used to store the multiplier m

            for (int iter = iter_y+1; iter < u_mtx._dim_y; iter ++)
            {
                u_mtx(sortedRowIndex[iter_x], iter) -= m * u_mtx(sortedRowIndex[iter_y], iter);
            }
        }
    }
}

template<class T> 
Matrix<double> Matrix<T>::GaussBackward(Matrix<double> b_mtx, Matrix<double> u_mtx, std::vector<int> sortedRowIndex)
{
// Solve Amat * Xmat = Bmat,
//  where Amat (n*n), Xmat (n*m), Bmat (n*m).
// In GaussForward, Amat = Lmat * u_mtx, 
//  where Lmat (n*n) and u_mtx (n*n),
//  we are trying to solve for u_mtx (upper triangle matrix).
// The actions (Lmat^-1) are stored in 
//  the lower triangle part of u_mtx (multipliers),
//  and the sortedRowIndex.
// In GaussBackward, we are trying to solve for Xmat,
//  where Xmat = u_mtx^(-1) * Lmat^(-1) * Bmat

    Matrix<double> res(b_mtx);
    double z; // right-hand-side variable 

    for (int iter = 0; iter < b_mtx._dim_y; iter ++)
    {
        for (int iter_y = 0; iter_y < u_mtx._dim_y-1; iter_y ++)
        {
            for (int iter_x = iter_y+1; iter_x < u_mtx._dim_x; iter_x ++)
            {
                b_mtx(sortedRowIndex[iter_x], iter) -= (u_mtx(sortedRowIndex[iter_x], iter_y) 
                                                      * b_mtx(sortedRowIndex[iter_y], iter));
            }
        }

        // sol to the last variable
        res(b_mtx._dim_x-1, iter) = b_mtx(sortedRowIndex[b_mtx._dim_x-1], iter) 
                                  / u_mtx(sortedRowIndex[u_mtx._dim_x-1], u_mtx._dim_y-1); 
         
        for (int iter_x = u_mtx._dim_x-2; iter_x > -1; iter_x --)
        {
            z = b_mtx(sortedRowIndex[iter_x], iter);
            for (int iter_y = iter_x+1; iter_y < u_mtx._dim_y; iter_y ++)
            {
                z -= (u_mtx(sortedRowIndex[iter_x], iter_y) * res(iter_y, iter));
            }
            res(iter_x, iter) = z / u_mtx(sortedRowIndex[iter_x], iter_x);
        }
        
    }

    return res;
}



#endif
