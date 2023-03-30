#include "myIO.h"

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

void WriteTracerPos (std::string file_name, std::vector<TracerInfo> const& object_list)
{
    std::cout << "Start writing!" << std::endl;

    std::ofstream outfile(file_name, std::ios::out);

    int n = object_list.size();
    Matrix<double> pos(3,1);

    if (n==0)
    {
        std::cout << "myIO::WriteObjectPos: n=0" << std::endl;
    }
    else 
    {
        for (int i = 0; i < n; i ++)
        {
            pos = object_list[i].GetCenterPos();
            outfile << pos(0,0) << ","
                    << pos(1,0) << "," 
                    << pos(2,0) << "\n";
        }
    }

    outfile.close();
    std::cout << "Finish writing!" << std::endl;
}

// .csv file 
// not recommended for large matrix
template<class T>
void LoadMatrix (std::string file_name, std::vector<std::vector<T>>& mtx)
{
    std::cout << "Start loading!" << std::endl;

    std::ifstream infile;
    infile.open(file_name);
    std::string line;
    T value;

    while (std::getline(infile, line, ','))
    {
        std::istringstream istream;
        istream.str(line);

        std::vector<T> row;
        while(istream >> value)
        {
            row.push_back(value);
        }
        mtx.push_back(row);
    }

    infile.close();
    std::cout << "Finish loading!" << std::endl;
};

}