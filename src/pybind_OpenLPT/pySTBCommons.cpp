#include <pybind11/pybind11.h>
#include "STBCommons.h"

namespace py = pybind11;
using namespace pybind11::literals;

void init_STBCommons(py::module &m)
{
    py::class_<PixelRange>(m, "PixelRange")
        .def(py::init<>())
        .def_readwrite("row_min", &PixelRange::row_min)
        .def_readwrite("row_max", &PixelRange::row_max)
        .def_readwrite("col_min", &PixelRange::col_min)
        .def_readwrite("col_max", &PixelRange::col_max)
        .def("setRowRange", &PixelRange::setRowRange)
        .def("setColRange", &PixelRange::setColRange)
        .def("setRange", &PixelRange::setRange)
        .def("getNumOfRow", &PixelRange::getNumOfRow)
        .def("getNumOfCol", &PixelRange::getNumOfCol)
        .def("to_dict", [](PixelRange const& self){
            return py::dict(
                "row_min"_a=self.row_min, 
                "row_max"_a=self.row_max, 
                "col_min"_a=self.col_min, 
                "col_max"_a=self.col_max
            );
        })
        .doc() = "PixelRange struct";

    py::class_<AxisLimit>(m, "AxisLimit")
        .def(py::init<>())
        .def(py::init<double, double, double, double, double, double>())
        .def_readwrite("x_min", &AxisLimit::x_min)
        .def_readwrite("x_max", &AxisLimit::x_max)
        .def_readwrite("y_min", &AxisLimit::y_min)
        .def_readwrite("y_max", &AxisLimit::y_max)
        .def_readwrite("z_min", &AxisLimit::z_min)
        .def_readwrite("z_max", &AxisLimit::z_max)
        .def("check", &AxisLimit::check)
        .def("to_dict", [](AxisLimit const& self){
            return py::dict(
                "x_min"_a=self.x_min, 
                "x_max"_a=self.x_max, 
                "y_min"_a=self.y_min, 
                "y_max"_a=self.y_max,
                "z_min"_a=self.z_min,
                "z_max"_a=self.z_max
            );
        })
        .doc() = "AxisLimit struct";

    py::enum_<ErrorTypeID>(m, "ErrorTypeID")
        .value("error_size", ErrorTypeID::error_size)
        .value("error_type", ErrorTypeID::error_type)
        .value("error_range", ErrorTypeID::error_range)
        .value("error_space", ErrorTypeID::error_space)
        .value("error_io", ErrorTypeID::error_io)
        .value("error_div0", ErrorTypeID::error_div0)
        .value("error_parallel", ErrorTypeID::error_parallel)
        .export_values();

    py::enum_<ObjectTypeID>(m, "ObjectTypeID")
        .value("type_tracer", ObjectTypeID::type_tracer)
        .value("type_bubble", ObjectTypeID::type_bubble)
        .value("type_filament", ObjectTypeID::type_filament)
        .export_values();

    py::enum_<FrameTypeID>(m, "FrameTypeID")
        .value("PREV_FRAME", FrameTypeID::PREV_FRAME)
        .value("CURR_FRAME", FrameTypeID::CURR_FRAME)
        .export_values();
}

