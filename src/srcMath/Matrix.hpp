#ifndef MATRIX_HPP
#define MATRIX_HPP

#include "Matrix.h"

// test 
template<class T>
Matrix<T>::Matrix (int dim_x, int dim_y, std::initializer_list<std::initializer_list<T>> mtx)
{
    create(dim_x, dim_y);

    for (int i = 0; i < _dim_x; i ++)
    {
        for (int j = 0; j < _dim_y; j ++)
        {
            _mtx[mapID(i,j)] = mtx.begin()[i].begin()[j];
        }
    }
}


//
// Constructor
//
template<class T> 
Matrix<T>::Matrix(int dim_x, int dim_y, T val)
{
    create(dim_x, dim_y);

    for (int i = 0; i < _n; i ++)
    {
        _mtx[i] = val;
    }
}

template<class T> 
Matrix<T>::Matrix (const Matrix<T>& mtx) 
{
    create(mtx._dim_x, mtx._dim_y);

    for (int i = 0; i < _n; i ++)
    {
        _mtx[i] = mtx._mtx[i];
    }
}

template<class T>
Matrix<T>::Matrix (std::string file_name)
{
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
        create(dim_x, dim_y);
        
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
        create(1, 1);
    }
}

// Deconstructor 
template<class T>
Matrix<T>::~Matrix ()
{
    clear();
}


// Create/clear space 
template<class T>
void Matrix<T>::clear()
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
void Matrix<T>::create(int dim_x, int dim_y)
{
    if (!_is_space)
    {
        _is_space = 1;

        _dim_x = dim_x;
        _dim_y = dim_y;
        _n = _dim_x * _dim_y;

        _mtx = new T [_n];
    }
    else
    {
        std::cerr << "The space for _mtx is not cleared: "
                  << "Cannot claim more space!" << std::endl;
        throw error_space;
    }
}


// Get id
template<class T>
inline int Matrix<T>::mapID(int id_x, int id_y) const
{
    return (id_x * _dim_y + id_y);
}


// Type transform
template<class T> 
Matrix<double> Matrix<T>::typeToDouble ()
{
    Matrix<double> res(_dim_x,_dim_y);
    for (int i = 0; i < _n; i ++)
    {
        res._mtx[i] = double(_mtx[i]);
    }
    return res;
}

template<class T> 
Matrix<int> Matrix<T>::typeToInt ()
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
    return _mtx[mapID(i,j)];
}

template<class T> 
T& Matrix<T>::operator() (int i, int j)
{
    return _mtx[mapID(i,j)];
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

//
// Get matrix info
//
template<class T> 
inline int Matrix<T>::getDimX() const
{
    return _dim_x;
}

template<class T> 
inline int Matrix<T>::getDimY() const
{
    return _dim_y;
}

template<class T> 
void Matrix<T>::print()
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
double Matrix<T>::norm()
{
    double res = 0.0;
    for (int i = 0; i < _n; i ++)
    {
        res += (_mtx[i] * _mtx[i]);
    }
    return std::sqrt(res);
}


//
// Matrix output 
//
template<class T>
void Matrix<T>::write (std::string file_name)
{
    std::cout << "\nStart writing!" << std::endl;

    std::ofstream outfile(file_name, std::ios::out);
    outfile.precision(8);

    for (int i = 0; i < _dim_x; i ++)
    {
        for (int j = 0; j < _dim_y-1; j ++)
        {
            outfile << _mtx[mapID(i,j)] << ",";
        }
        outfile << _mtx[mapID(i,_dim_y-1)] << "\n";
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
        clear();
        create(mtx._dim_x, mtx._dim_y);
    }

    for (int i = 0; i < _n; i ++)
    {
        _mtx[i]= mtx._mtx[i];
    }

    return *this;
}

template<class T>
bool Matrix<T>::operator== (Matrix<T> const& mtx1)
{
    if (_dim_x != mtx1._dim_x || _dim_y != mtx1._dim_y)
    {
        return false;
    }

    for (int i = 0; i < _n; i ++)
    {
        if (std::fabs(_mtx[i] - mtx1._mtx[i]) > SMALLNUMBER)
        {
            return false;
        }
    }

    return true;
}

template<class T>
bool Matrix<T>::operator!= (Matrix<T> const& mtx1)
{
    if (_dim_x != mtx1._dim_x || _dim_y != mtx1._dim_y)
    {
        return true;
    }

    for (int i = 0; i < _n; i ++)
    {
        if (std::fabs(_mtx[i] - mtx1._mtx[i]) > SMALLNUMBER)
        {
            return true;
        }
    }

    return false;
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
                res(iter_x,iter_y) += (_mtx[mapID(iter_x,iter)] * mtx(iter,iter_y));
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
                res(iter_x,iter_y) += (_mtx[mapID(iter_x,iter)] * mtx(iter,iter_y));
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


// Get Row 
template<class T>
std::vector<T> Matrix<T>::getRow(int row_i) const
{
    std::vector<T> res(_dim_y);
    for (int j = 0; j < _dim_y; j ++)
    {
        res[j] = _mtx[mapID(row_i,j)];
    }
    return res;
}

// Get Col 
template<class T>
std::vector<T> Matrix<T>::getCol(int col_j) const
{
    std::vector<T> res(_dim_x);
    for (int i = 0; i < _dim_x; i ++)
    {
        res[i] = _mtx[mapID(i,col_j)];
    }
    return res;
}

//
// Matrix manipulation
//
template<class T>
Matrix<T> Matrix<T>::transpose ()
{
    Matrix<T> res(_dim_y, _dim_x, 0);
    for (int i = 0; i < _dim_x; i ++)
    {
        for (int j = 0; j < _dim_y; j ++)
        {
            res(j,i) = _mtx[mapID(i,j)];
        }
    }
    return res;
}


#endif
