#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include "Matrix.h"

namespace py = pybind11;

template <class T>
py::array_t<T> matrix_to_numpy(Matrix<T> const& mat) 
{
    // Define the shape of the array (rows, cols)
    size_t n_row = mat.getDimRow();
    size_t n_col = mat.getDimCol();
    std::vector<size_t> shape = {n_row, n_col};

    // Create a NumPy array with the same size and data type as the matrix
    py::array_t<T> array(shape);

    // Copy data from matrix to NumPy array
    auto buf = array.request();
    T* ptr = static_cast<T*>(buf.ptr);
    std::copy(mat.data(), mat.data() + (n_row*n_col), ptr);

    return array;
}

void init_Matrix(py::module &m) 
{
    py::class_<Matrix<double>>(m, "Matrix")
        .def(py::init<int, int, double>())
        .def("operator=", &Matrix<double>::operator=)
        .def("__setitem__", [](Matrix<double> &self, int id, float val) {self[id] = val;})
        .def("__setitem__", [](Matrix<double> &self, std::pair<int,int> index, float val) {
            self(index.first,index.second) = val;
        })
        .def("__getitem__", [](Matrix<double> const&self, int id) {return self[id];})
        .def("__getitem__", [](Matrix<double> const&self, std::pair<int,int> index) {
            return self(index.first,index.second);
        });
    
    m.def("matrix_to_numpy", &matrix_to_numpy<double>, "Convert a Matrix<double> to a NumPy array");
}

