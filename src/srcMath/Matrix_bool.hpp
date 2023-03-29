#ifndef MATRIXBOOL_HPP
#define MATRIXBOOL_HPP

#include "Matrix_bool.h"

//                           
// Constructor/Deconstructor 
//                           
Matrix_bool::Matrix_bool (int dim_x, int dim_y)
{
    CreateMtx(dim_x, dim_y);
}

Matrix_bool::Matrix_bool (std::string file_name)
{
    std::cout << "Start loading!" << std::endl;

    std::string line;
    bool value;

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
        CreateMtx(dim_x, dim_y);
        
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
        CreateMtx(1,1);
    }

    std::cout << "Finish loading!" << std::endl;
}

Matrix_bool::~Matrix_bool ()
{
    ClearMtx();
}


//
// Create/Clear space
//
void Matrix_bool::CreateMtx(int dim_x, int dim_y)
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
                _mtx = new bool* [_dim_x];

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
                    _mtx[i] = new bool [_dim_y];

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
                _mtx[i][j] = false;
            }
        }
    }
    else
    {
        std::cerr << "The space for _mtx is not cleared: "
                  << "Cannot claim more space!" << std::endl;
        throw error_space;
    }
}

void Matrix_bool::ClearMtx()
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


//
// Get matrix info
//
int Matrix_bool::GetDimX() const
{
    return _dim_x;
}

int Matrix_bool::GetDimY() const
{
    return _dim_y;
}

void Matrix_bool::Print()
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


//
// Get/Assign value
//
bool Matrix_bool::operator() (int i, int j) const
{
    return _mtx[i][j];
}
 
bool& Matrix_bool::operator() (int i, int j)
{
    return _mtx[i][j];
}

bool Matrix_bool::operator[] (int i) const 
{
    return _mtx[i][0];
}

bool& Matrix_bool::operator[] (int i)
{
    return _mtx[i][0];
}

std::vector<bool> Matrix_bool::GetRow(int i)
{
    std::vector<bool> row(_dim_y, 0);
    for (int iter = 0; iter < _dim_y; iter ++)
    {
        row[iter] = _mtx[i][iter];
    }
    return row;
}

std::vector<bool> Matrix_bool::GetCol(int j)
{
    std::vector<bool> col(_dim_x, 0);
    for (int iter = 0; iter < _dim_x; iter ++)
    {
        col[iter] = _mtx[iter][j];
    }
    return col;
}


//
// Matrix output 
//
void Matrix_bool::WriteMatrix (std::string file_name)
{
    std::cout << "Start writing!" << std::endl;

    std::ofstream outfile(file_name, std::ios::out);

    outfile << "Matrix: " << "(row,col)=" 
            << "(" << _dim_x << "," 
            << _dim_y << ")" << "\n";

    bool value;
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


#endif