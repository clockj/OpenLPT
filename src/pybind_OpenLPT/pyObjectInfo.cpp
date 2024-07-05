#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <iostream>
#include <string>
#include <sstream>
#include <fstream>

#include "ObjectInfo.h"

namespace py = pybind11;


void init_ObjectInfo(py::module &m) 
{
    py::class_<Object2D>(m, "Object2D")
        .def(py::init<>())
        .def(py::init<Object2D const&>())
        .def(py::init<Pt2D const&>())
        .def_readwrite("_pt_center", &Object2D::_pt_center)
        .doc() = "Object2D class";

    py::class_<Tracer2D, Object2D>(m, "Tracer2D")
        .def(py::init<>())
        .def(py::init<Tracer2D const&>())
        .def(py::init<Pt2D const&>())
        .def_readwrite("_r_px", &Tracer2D::_r_px)
        .doc() = "Tracer2D class";

    py::class_<Object3D>(m, "Object3D")
        .def(py::init<>())
        .def(py::init<Object3D const&>())
        .def(py::init<Pt3D const&>())
        .def_readwrite("_pt_center", &Object3D::_pt_center)
        .def_readwrite("_is_tracked", &Object3D::_is_tracked)
        .doc() = "Object3D class";

    py::class_<Tracer3D, Object3D>(m, "Tracer3D")
        .def(py::init<>())
        .def(py::init<Tracer3D const&>())
        .def(py::init<Pt3D const&>())
        .def_readwrite("_n_2d", &Tracer3D::_n_2d)
        .def_readwrite("_error", &Tracer3D::_error)
        .def_readwrite("_camid_list", &Tracer3D::_camid_list)
        .def_readwrite("_tr2d_list", &Tracer3D::_tr2d_list)
        .def("addTracer2D", [](Tracer3D& self, Tracer2D const& tracer2d, int cam_id){
            self.addTracer2D(tracer2d, cam_id);
        })
        .def("addTracer2D", [](Tracer3D& self, std::vector<Tracer2D> const& tracer2d_list, std::vector<int> const& camid_list){
            self.addTracer2D(tracer2d_list, camid_list);
        })
        .def("removeTracer2D", [](Tracer3D& self, int cam_id){
            self.removeTracer2D(cam_id);
        })
        .def("removeTracer2D", [](Tracer3D& self, std::vector<int> const& camid_list){
            self.removeTracer2D(camid_list);
        })
        .def("clearTracer2D", &Tracer3D::clearTracer2D)
        .def("updateTracer2D", [](Tracer3D& self, Tracer2D const& tracer2d, int cam_id){
            self.updateTracer2D(tracer2d, cam_id);
        })
        .def("updateTracer2D", [](Tracer3D& self, std::vector<Tracer2D> const& tracer2d_list, std::vector<int> const& camid_list){
            self.updateTracer2D(tracer2d_list, camid_list);
        })
        .def("projectObject2D", &Tracer3D::projectObject2D)
        .def("getTracer2D", [](Tracer3D& self, int cam_id){
            Tracer2D tracer2d;
            self.getTracer2D(tracer2d, cam_id);
            return tracer2d;
        })
        .def("saveObject3D", [](Tracer3D& self, std::string& file, int n_cam_all, bool is_append=true){
            std::ofstream output(file, is_append ? std::ios::app : std::ios::out);
            self.saveObject3D(output, n_cam_all);
        }, py::arg("file"), py::arg("n_cam_all"), py::arg("is_append")=true)
        .doc() = "Tracer3D class";
}