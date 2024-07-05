#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "ImageIO.h"

namespace py = pybind11;


void init_ImageIO(py::module &m) 
{
    py::class_<ImageParam>(m, "ImageParam")
        .def(py::init<>())
        .def_readwrite("n_row", &ImageParam::n_row)
        .def_readwrite("n_col", &ImageParam::n_col)
        .def_readwrite("bits_per_sample", &ImageParam::bits_per_sample)
        .def_readwrite("n_channel", &ImageParam::n_channel)
        .doc() = "ImageParam struct";

    py::class_<ImageIO>(m, "ImageIO")
        .def(py::init<>())
        .def(py::init<ImageIO const&>())
        .def("loadImgPath", &ImageIO::loadImgPath)
        .def("loadImg", &ImageIO::loadImg)
        .def("saveImg", &ImageIO::saveImg)
        .def("setImgParam", &ImageIO::setImgParam)
        .def("getImgParam", &ImageIO::getImgParam)
        .def("getCurrImgID", &ImageIO::getCurrImgID)
        .doc() = "ImageIO class";
}