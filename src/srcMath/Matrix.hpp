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
    : Matrix(mtx.size(), mtx[0].size(), 0)
{
    for (int i = 0; i < _dim_x; i ++)
    {
        if (mtx[i].size() != _dim_y)
        {
            std::cerr << "Matrix::Matrix: " 
                      << "_dim_y=" << _dim_y << ","
                      << "mtx[" << i << "].size()=" << mtx[i].size()
                      << std::endl;
            throw error_size;
        }

        for (int j = 0; j < _dim_y; j ++)
        {
            _mtx[i][j] = mtx[i][j];
        }
    }
}

template<class T> 
Matrix<T>::Matrix (const Matrix<T>& mtx) 
    : Matrix(mtx._dim_x, mtx._dim_y, 0)
{
    for (int i = 0; i < _dim_x; i ++)
    {
        for (int j = 0; j < _dim_y; j ++)
        {
            _mtx[i][j] = mtx._mtx[i][j];
        }
    }
}

template<class T>
Matrix<T>::Matrix (const Matrix<T>& vec_1, const Matrix<T>& vec_2)
    : Matrix(vec_1._dim_x, vec_1._dim_y, 0)
{
    if (vec_2._dim_x != _dim_x || vec_2._dim_y != _dim_y)
    {
        std::cerr << "Matrix::Matrix: the sizes of vec_1 and vec_2 are different" 
                  << std::endl;
        throw error_size;
    }

    T mag = 0;
    T value = 0;
    for (int i = 0; i < _dim_x; i ++)
    {
        for (int j = 0; j < _dim_y; j ++)
        {
            value = vec_1._mtx[i][j] - vec_2._mtx[i][j];

            _mtx[i][j] = value;
            mag += value * value;
        }
    }

    if (mag < MAGSMALLNUMBER)
    {
        std::cerr << "Matrix::Matrix: vec_1 and vec_2 are too close"
                  << std::endl;
        throw error_div0;
    }

    mag = T (std::sqrt(mag));
    for (int i = 0; i < _dim_x; i ++)
    {
        for (int j = 0; j < _dim_y; j ++)
        {
            _mtx[i][j] /= mag;
        }
    }
}

template<class T> 
Matrix<T>::Matrix (T m11, T m12,
                   T m21, T m22)
    : Matrix(2, 2, 0)
{
    _mtx[0][0] = m11; _mtx[0][1] = m12;
    _mtx[1][0] = m21; _mtx[1][1] = m22;
}

template<class T> 
Matrix<T>::Matrix (T m11, T m12, T m13,
                   T m21, T m22, T m23,
                   T m31, T m32, T m33)
    : Matrix(3, 3, 0)
{    
    _mtx[0][0] = m11; _mtx[0][1] = m12; _mtx[0][2] = m13;
    _mtx[1][0] = m21; _mtx[1][1] = m22; _mtx[1][2] = m23;
    _mtx[2][0] = m31; _mtx[2][1] = m32; _mtx[2][2] = m33;
}

template<class T>
Matrix<T>::Matrix (int dim_x, int dim_y, T v1, T v2, T v3)
    : Matrix(3,1)
{
    if (dim_x != 3 || dim_y != 1)
    {
        std::cerr << "The dimension does not satisfied for "
                  << "creating a 3*1 vector!"
                  << std::endl;
        throw error_size;
    }

    _mtx[0][0] = v1;
    _mtx[1][0] = v2;
    _mtx[2][0] = v3;
}

template<class T> 
Matrix<T>::Matrix (int dim)
    : Matrix(dim, dim, 0)
{
    for (int i = 0; i < _dim_x; i ++)
    {
        _mtx[i][i] = 1;
    }
}

template<class T>
Matrix<T>::Matrix (std::string file_name)
{
    std::cout << "Start loading!" << std::endl;

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

            int j = 0;
            while (istream_data >> value)
            {
                _mtx[i][j] = value;
                j ++;
                if (istream_data.peek() == ',')
                {
                    istream_data.ignore();
                }
            }
            if (j < _dim_y)
            {
                for (int z = j; z < _dim_y; z ++)
                {
                    _mtx[i][z] = 0;
                }
            }
            
        }

        infile_data.close();
    }
    else 
    {
        CreateMtx(1,1,0);
    }

    std::cout << "Finish loading!" << std::endl;
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
        for (int i = 0; i < _dim_x; i ++)
        {
            delete[] _mtx[i];
        }
        delete[] _mtx;
        _is_space = 0;
    }
}

template<class T>
void Matrix<T>::ReCreateCol (int dim_y)
{
    if (_is_space && dim_y > 0)
    {
        _dim_y = dim_y;
        for (int i = 0; i < _dim_x; i ++)
        {
            delete[] _mtx[i];
            _mtx[i] = new T [_dim_y];

            for (int j = 0; j < _dim_y; j ++)
            {
                _mtx[i][j] = 0;
            }
        }
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

        int temp = 0;
        while (temp < 5)
        {
            try 
            {
                _mtx = new T* [_dim_x];

                break;
            }
            catch (...)
            {    
                if (temp == 4)
                {
                    std::cout << "CreateMtx: dim_x no enough space\n";
                    throw;
                }
            }
            temp += 1;
        }
        for (int i = 0; i < _dim_x; i ++)
        {
            int temp = 0;
            while (temp < 5)
            {
                try 
                {
                    _mtx[i] = new T [_dim_y];

                    break;
                }
                catch (...)
                {    
                    if (temp == 4)
                    {
                        std::cout << "CreateMtx: dim_y no enough space\n";
                        throw;
                    }
                }
                temp += 1;
            }
            
            for (int j = 0; j < _dim_y; j ++)
            {
                _mtx[i][j] = val;
            }
        }

        // _mtx = new T* [_dim_x];
        // for (int i = 0; i < _dim_x; i ++)
        // {
        //     _mtx[i] = new T [_dim_y];
        //     for (int j = 0; j < _dim_y; j ++)
        //     {
        //         _mtx[i][j] = val;
        //     }
        // }
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


// Type transform
template<class T> 
Matrix<double> Matrix<T>::TypeToDouble ()
{
    Matrix<double> res(_dim_x, _dim_y);
    for (int iter_x = 0; iter_x < _dim_x; iter_x ++)
    {
        for (int iter_y = 0; iter_y < _dim_y; iter_y ++)
        {
            res(iter_x, iter_y) = double(_mtx[iter_x][iter_y]);
        }
    }
    return res;
}

template<class T> 
Matrix<int> Matrix<T>::TypeToInt ()
{
    Matrix<int> res(_dim_x, _dim_y);
    for (int iter_x = 0; iter_x < _dim_x; iter_x ++)
    {
        for (int iter_y = 0; iter_y < _dim_y; iter_y ++)
        {
            res(iter_x, iter_y) = int(_mtx[iter_x][iter_y]);
        }
    }
    return res;
}



//
// Get/Assign value
//
template<class T> 
T Matrix<T>::operator() (int i, int j) const
{
    return _mtx[i][j];
}

template<class T> 
T& Matrix<T>::operator() (int i, int j)
{
    return _mtx[i][j];
}

template<class T>
T Matrix<T>::operator[] (int i) const 
{
    return _mtx[i][0];
}

template<class T>
T& Matrix<T>::operator[] (int i)
{
    return _mtx[i][0];
}

template<class T> 
std::vector<T> Matrix<T>::GetRow(int i)
{
    std::vector<T> row(_dim_y, 0);
    for (int iter = 0; iter < _dim_y; iter ++)
    {
        row[iter] = _mtx[i][iter];
    }
    return row;
}

template<class T> 
std::vector<T> Matrix<T>::GetCol(int j)
{
    std::vector<T> col(_dim_x, 0);
    for (int iter = 0; iter < _dim_x; iter ++)
    {
        col[iter] = _mtx[iter][j];
    }
    return col;
}


//
// Get matrix info
//
template<class T> 
int Matrix<T>::GetDimX() const
{
    return _dim_x;
}

template<class T> 
int Matrix<T>::GetDimY() const
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
    
    for (int iter_x = 0; iter_x < _dim_x; iter_x ++)
    {
        std::cout << "| ";
        for (int iter_y = 0; iter_y < _dim_y; iter_y ++)
        {
            std::cout << std::scientific 
                      << std::setprecision(2)  
                      << std::setfill(' ') 
                      << std::left 
                      << _mtx[iter_x][iter_y] << " ";
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
    for (int iter_x = 0; iter_x < _dim_x; iter_x ++)
    {
        for (int iter_y = 0; iter_y < _dim_y; iter_y ++)
        {
            res += _mtx[iter_x][iter_y] * _mtx[iter_x][iter_y];
        }
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
    for (int i = 0; i < _dim_x; i ++)
    {
        for (int j = 0; j < _dim_y; j ++)
        {
            res += _mtx[i][j] * mtx._mtx[i][j];
        }
    }

    return res;
}


//
// Matrix output 
//
template<class T>
void Matrix<T>::WriteMatrix (std::string file_name)
{
    std::cout << "Start writing!" << std::endl;

    std::ofstream outfile(file_name, std::ios::out);

    outfile << "Matrix: " << "(row,col)=" 
            << "(" << _dim_x << "," 
            << _dim_y << ")" << "\n";

    T value;
    for (int i = 0; i < _dim_x; i ++)
    {
        for (int j = 0; j < _dim_y-1; j ++)
        {
            value = _mtx[i][j];
            outfile << value << ",";
        }
        outfile << _mtx[i][_dim_y-1] << "\n";
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
        CreateMtx(mtx._dim_x, mtx._dim_y, 0);
    }

    for (int i = 0; i < _dim_x; i ++)
    {
        for (int j = 0; j < _dim_y; j ++)
        {
            _mtx[i][j] = mtx._mtx[i][j];
        }
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
    for (int iter_x = 0; iter_x < _dim_x; iter_x ++)
    {
        for (int iter_y = 0; iter_y < _dim_y; iter_y ++)
        {
            res._mtx[iter_x][iter_y] = _mtx[iter_x][iter_y] + mtx._mtx[iter_x][iter_y];
        }
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
    for (int iter_x = 0; iter_x < _dim_x; iter_x ++)
    {
        for (int iter_y = 0; iter_y < _dim_y; iter_y ++)
        {
            _mtx[iter_x][iter_y] += mtx._mtx[iter_x][iter_y];
        }
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

    for (int iter_x = 0; iter_x < _dim_x; iter_x ++)
    {
        for (int iter_y = 0; iter_y < _dim_y; iter_y ++)
        {
            res._mtx[iter_x][iter_y] = _mtx[iter_x][iter_y] - mtx._mtx[iter_x][iter_y];
        }
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
    for (int iter_x = 0; iter_x < _dim_x; iter_x ++)
    {
        for (int iter_y = 0; iter_y < _dim_y; iter_y ++)
        {
            _mtx[iter_x][iter_y] -= mtx._mtx[iter_x][iter_y];
        }
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
                res._mtx[iter_x][iter_y] += (_mtx[iter_x][iter] * mtx._mtx[iter][iter_y]);
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
                res._mtx[iter_x][iter_y] += (_mtx[iter_x][iter] * mtx._mtx[iter][iter_y]);
            }
        }
    }
    
    if (_dim_y != res._dim_y)
    {
        ReCreateCol(res._dim_y);
    }

    for (int i = 0; i < _dim_x; i ++)
    {
        for (int j = 0; j < _dim_y; j ++)
        {
            _mtx[i][j] = res._mtx[i][j];
        }
    }

    return *this;
}


template<class T> 
Matrix<T> Matrix<T>::operator* (T ratio)
{
    Matrix<T> res(*this);
    for (int iter_x = 0; iter_x < _dim_x; iter_x ++)
    {
        for (int iter_y = 0; iter_y < _dim_y; iter_y ++)
        {
            res._mtx[iter_x][iter_y] = _mtx[iter_x][iter_y] * ratio;
        }
    }

    return res;
}

template<class T> 
Matrix<T>& Matrix<T>::operator*= (T ratio)
{
    for (int iter_x = 0; iter_x < _dim_x; iter_x ++)
    {
        for (int iter_y = 0; iter_y < _dim_y; iter_y ++)
        {
            _mtx[iter_x][iter_y] *= ratio;   
        }
    }

    return *this;
}

template<class T> 
Matrix<T> Matrix<T>::operator/ (T ratio)
{
    Matrix<T> res(*this);
    for (int iter_x = 0; iter_x < _dim_x; iter_x ++)
    {
        for (int iter_y = 0; iter_y < _dim_y; iter_y ++)
        {
            res._mtx[iter_x][iter_y] = _mtx[iter_x][iter_y] / ratio;
        }
    }

    return res;
}

template<class T> 
Matrix<T>& Matrix<T>::operator/= (T ratio)
{
    for (int iter_x = 0; iter_x < _dim_x; iter_x ++)
    {
        for (int iter_y = 0; iter_y < _dim_y; iter_y ++)
        {
            _mtx[iter_x][iter_y] /= ratio;   
        }
    }

    return *this;
}


template<class T> 
void Matrix<T>::PiecewiseMultiply (Matrix<T> const& mtx)
{
    for (int iter_x = 0; iter_x < _dim_x; iter_x ++)
    {
        for (int iter_y = 0; iter_y < _dim_y; iter_y ++)
        {
            _mtx[iter_x][iter_y] *= mtx._mtx[iter_x][iter_y];
        }
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
            res(iter_y, iter_x) = _mtx[iter_x][iter_y];
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
        res(iter_x, iter_x) = _mtx[iter_x][iter_x];
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
        res += _mtx[iter_x][iter_x];
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

    double det =   _mtx[0][0] * ( _mtx[1][1]*_mtx[2][2] - _mtx[1][2]*_mtx[2][1] )
                 - _mtx[0][1] * ( _mtx[1][0]*_mtx[2][2] - _mtx[1][2]*_mtx[2][0] ) 
                 + _mtx[0][2] * ( _mtx[1][0]*_mtx[2][1] - _mtx[1][1]*_mtx[2][0] );
    res._mtx[0][0] =   (_mtx[1][1]*_mtx[2][2] - _mtx[1][2]*_mtx[2][1]) / det;
    res._mtx[0][1] = - (_mtx[1][0]*_mtx[2][2] - _mtx[1][2]*_mtx[2][0]) / det;
    res._mtx[0][2] =   (_mtx[1][0]*_mtx[2][1] - _mtx[1][1]*_mtx[2][0]) / det;

    res._mtx[1][0] = - (_mtx[0][1]*_mtx[2][2] - _mtx[0][2]*_mtx[2][1]) / det;
    res._mtx[1][1] =   (_mtx[0][0]*_mtx[2][2] - _mtx[0][2]*_mtx[2][0]) / det;
    res._mtx[1][2] = - (_mtx[0][0]*_mtx[2][1] - _mtx[0][1]*_mtx[2][0]) / det;

    res._mtx[2][0] =   (_mtx[0][1]*_mtx[1][2] - _mtx[0][2]*_mtx[1][1]) / det;
    res._mtx[2][1] = - (_mtx[0][0]*_mtx[1][2] - _mtx[0][2]*_mtx[1][0]) / det;
    res._mtx[2][2] =   (_mtx[0][0]*_mtx[1][1] - _mtx[0][1]*_mtx[1][0]) / det;

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
                b_mtx(sortedRowIndex[iter_x], iter) -= u_mtx(sortedRowIndex[iter_x], iter_y) 
                                                     * b_mtx(sortedRowIndex[iter_y], iter);
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
                z -= u_mtx(sortedRowIndex[iter_x], iter_y) * res(iter_y, iter);
            }
            res(iter_x, iter) = z / u_mtx(sortedRowIndex[iter_x], iter_x);
        }
        
    }

    return res;
}



#endif
