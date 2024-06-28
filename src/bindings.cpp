#include <pybind11/pybind11.h>

namespace py = pybind11;

void init_Matrix(py::module &);

PYBIND11_MODULE(pyOpenLPT, m) 
{
    py::module m_math = m.def_submodule("math", "Math module");

    init_Matrix(m_math);
}