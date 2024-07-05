#include <pybind11/pybind11.h>
#include "STBCommons.h"

namespace py = pybind11;


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
        .doc() = "PixelRange struct";

    py::class_<AxisLimit>(m, "AxisLimit")
        .def(py::init<>())
        .def_readwrite("x_min", &AxisLimit::x_min)
        .def_readwrite("x_max", &AxisLimit::x_max)
        .def_readwrite("y_min", &AxisLimit::y_min)
        .def_readwrite("y_max", &AxisLimit::y_max)
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

