#ifndef MYIO_H
#define MYIO_H

#include <iostream>
#include <string>
#include <sstream>
#include <fstream>

#include "Matrix.h"
#include "ObjectInfo.h"

namespace myIO
{

template<class T>
void WriteMatrix (std::string file_name, std::vector<std::vector<T>> const& mtx)
{
    std::cout << "Start writing!" << std::endl;

    std::ofstream outfile(file_name, std::ios::out);

    int dim_x = mtx.size();

    if (dim_x==0)
    {
        std::cout << "myIO::WriteMatrix: dim_x=0" << std::endl;
    }
    else 
    {
        for (int i = 0; i < dim_x; i ++)
        {
            int dim_y = mtx[i].size();

            if (dim_y == 0)
            {
                outfile << "\n";
            }
            else 
            {
                for (int j = 0; j < dim_y-1; j ++)
                {
                    outfile << mtx[i][j] << ",";
                }
                outfile << mtx[i][dim_y-1] << "\n";
            }
        }
    }

    outfile.close();
    std::cout << "Finish writing!" << std::endl;
};

template<class T>
void WriteMatrix (std::string file_name, std::vector<T> const& mtx)
{
    std::cout << "Start writing!" << std::endl;

    std::ofstream outfile(file_name, std::ios::out);

    int dim_x = mtx.size();

    if (dim_x==0)
    {
        std::cout << "myIO::WriteMatrix: dim_x=0" << std::endl;
    }
    else 
    {
        for (int i = 0; i < dim_x-1; i ++)
        {
            outfile << mtx[i] << ",";
        }
        outfile << mtx[dim_x-1];
    }

    outfile.close();
    std::cout << "Finish writing!" << std::endl;
};

void WriteTracerPos (std::string file_name, std::vector<TracerInfo> const& object_list);

// .csv file 
// not recommended for large matrix
template<class T>
void LoadMatrix (std::string file_name, std::vector<std::vector<T>>& mtx);

}


#endif