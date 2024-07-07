#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <iostream>
#include <string>
#include <sstream>
#include <fstream>

#include "myMATH.h"
#include "Matrix.h"
#include "Camera.h"

namespace py = pybind11;
using namespace pybind11::literals;

void init_Camera(py::module &m) 
{
    py::class_<PinholeParam>(m, "PinholeParam")
        .def(py::init<>())
        .def_readwrite("n_row", &PinholeParam::n_row)
        .def_readwrite("n_col", &PinholeParam::n_col)
        .def_readwrite("cam_mtx", &PinholeParam::cam_mtx)
        .def_readwrite("is_distorted", &PinholeParam::is_distorted)
        .def_readwrite("n_dist_coeff", &PinholeParam::n_dist_coeff)
        .def_readwrite("dist_coeff", &PinholeParam::dist_coeff)
        .def_readwrite("r_mtx", &PinholeParam::r_mtx)
        .def_readwrite("t_vec", &PinholeParam::t_vec)
        .def_readwrite("r_mtx_inv", &PinholeParam::r_mtx_inv)
        .def_readwrite("t_vec_inv", &PinholeParam::t_vec_inv)
        .def("to_dict", [](PinholeParam const& self){
            return py::dict(
                "n_row"_a=self.n_row, 
                "n_col"_a=self.n_col, 
                "cam_mtx"_a=self.cam_mtx, 
                "is_distorted"_a=self.is_distorted, 
                "n_dist_coeff"_a=self.n_dist_coeff, 
                "dist_coeff"_a=self.dist_coeff, 
                "r_mtx"_a=self.r_mtx, 
                "t_vec"_a=self.t_vec, 
                "r_mtx_inv"_a=self.r_mtx_inv, 
                "t_vec_inv"_a=self.t_vec_inv
            );
        })
        .doc() = "PinholeParam struct";

    py::enum_<RefPlane>(m, "RefPlane")
        .value("REF_X", RefPlane::REF_X)
        .value("REF_Y", RefPlane::REF_Y)
        .value("REF_Z", RefPlane::REF_Z)
        .export_values();

    py::class_<PolyParam>(m, "PolyParam")
        .def(py::init<>())
        .def_readwrite("n_row", &PolyParam::n_row)
        .def_readwrite("n_col", &PolyParam::n_col)
        .def_readwrite("ref_plane", &PolyParam::ref_plane)
        .def_readwrite("plane", &PolyParam::plane)
        .def_readwrite("n_coeff", &PolyParam::n_coeff)
        .def_readwrite("u_coeffs", &PolyParam::u_coeffs)
        .def_readwrite("du_coeffs", &PolyParam::du_coeffs)
        .def_readwrite("v_coeffs", &PolyParam::v_coeffs)
        .def_readwrite("dv_coeffs", &PolyParam::dv_coeffs)
        .def("to_dict", [](PolyParam const& self){
            return py::dict(
                "n_row"_a=self.n_row, 
                "n_col"_a=self.n_col, 
                "ref_plane"_a=self.ref_plane, 
                "plane"_a=self.plane, 
                "n_coeff"_a=self.n_coeff, 
                "u_coeffs"_a=self.u_coeffs, 
                "du_coeffs"_a=self.du_coeffs, 
                "v_coeffs"_a=self.v_coeffs, 
                "dv_coeffs"_a=self.dv_coeffs
            );
        })
        .doc() = "PolyParam struct";

    py::enum_<CameraType>(m, "CameraType")
        .value("PINHOLE", CameraType::PINHOLE)
        .value("POLYNOMIAL", CameraType::POLYNOMIAL)
        .export_values();

    py::class_<Camera>(m, "Camera")
        .def_readwrite("_type", &Camera::_type)
        .def_readwrite("_pinhole_param", &Camera::_pinhole_param)
        .def_readwrite("_poly_param", &Camera::_poly_param)
        .def(py::init<>())
        .def(py::init<const Camera&>())
        .def(py::init<std::istream&>())
        .def(py::init<std::string>())
        .def("loadParameters", [](Camera &self, std::istream& is) {
            self.loadParameters(is);
        })
        .def("loadParameters", [](Camera &self, std::string filename) {
            self.loadParameters(filename);
        })
        .def("saveParameters", &Camera::saveParameters)
        .def("rmtxTorvec", &Camera::rmtxTorvec)
        .def("getNRow", &Camera::getNRow)
        .def("getNCol", &Camera::getNCol)
        .def("project", &Camera::project)
        .def("worldToUndistImg", &Camera::worldToUndistImg)
        .def("distort", &Camera::distort)
        .def("polyProject", &Camera::polyProject)
        .def("lineOfSight", &Camera::lineOfSight)
        .def("undistort", &Camera::undistort)
        .def("pinholeLine", &Camera::pinholeLine)
        .def("polyImgToWorld", &Camera::polyImgToWorld)
        .def("polyLineOfSight", &Camera::polyLineOfSight)
        .def("to_dict", [](Camera const& self){
            return py::dict(
                "_type"_a=self._type, 
                "_pinhole_param"_a=self._pinhole_param, 
                "_poly_param"_a=self._poly_param
            );
        })
        .doc() = "Camera class";

    py::class_<CamList>(m, "CamList")
        .def(py::init<>())
        .def(py::init([](std::vector<Camera> const& cam_list, std::vector<int> const& intensity_max, std::vector<int> const& useid_list){
            CamList cam_list_all;
            cam_list_all.cam_list = cam_list;
            cam_list_all.intensity_max = intensity_max;
            cam_list_all.useid_list = useid_list;
            return cam_list_all;
        }))
        .def_readwrite("cam_list", &CamList::cam_list)
        .def_readwrite("intensity_max", &CamList::intensity_max)
        .def_readwrite("useid_list", &CamList::useid_list)
        .def("to_dict", [](CamList const& self){
            return py::dict(
                "cam_list"_a=self.cam_list, 
                "intensity_max"_a=self.intensity_max, 
                "useid_list"_a=self.useid_list
            );
        })
        .doc() = "CamList struct";
}