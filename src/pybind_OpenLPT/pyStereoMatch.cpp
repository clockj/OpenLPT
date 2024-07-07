#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "StereoMatch.h"

namespace py = pybind11;
using namespace pybind11::literals;

void init_StereoMatch(py::module& m)
{
    py::class_<ObjIDMap>(m, "ObjIDMap")
        .def(py::init<>())
        .def("config", &ObjIDMap::config)
        .def("mapImgID", &ObjIDMap::mapImgID)
        .doc() = "ObjIDMap class";

    py::class_<SMParam>(m, "SMParam")
        .def(py::init<>())
        .def_readwrite("tor_2d", &SMParam::tor_2d)
        .def_readwrite("tor_3d", &SMParam::tor_3d)
        .def_readwrite("n_thread", &SMParam::n_thread)
        .def_readwrite("check_id", &SMParam::check_id)
        .def_readwrite("check_radius", &SMParam::check_radius)
        .def_readwrite("is_delete_ghost", &SMParam::is_delete_ghost)
        .def_readwrite("is_update_inner_var", &SMParam::is_update_inner_var)
        .def("to_dict", [](SMParam const& self){
            return py::dict(
                "tor_2d"_a=self.tor_2d, 
                "tor_3d"_a=self.tor_3d, 
                "n_thread"_a=self.n_thread, 
                "check_id"_a=self.check_id, "check_radius"_a=self.check_radius, 
                "is_delete_ghost"_a=self.is_delete_ghost, "is_update_inner_var"_a=self.is_update_inner_var
            );
        })
        .doc() = "SMParam struct";

    py::class_<StereoMatch>(m, "StereoMatch")
        .def(py::init<SMParam const&, CamList const&>())
        .def("clearAll", &StereoMatch::clearAll)
        .def("match", [](StereoMatch& self, std::vector<std::vector<Tracer2D>> const& obj2d_list){
            std::vector<Tracer3D> obj3d_list;
            self.match(obj3d_list, obj2d_list);
            return obj3d_list;
        })
        .def("saveObjInfo", [](StereoMatch& self, std::string path, std::vector<Tracer3D> const& obj3d_list){
            self.saveObjInfo(path, obj3d_list);
        })
        .def("saveObjIDMatchList", &StereoMatch::saveObjIDMatchList)
        .def_readwrite("_param", &StereoMatch::_param)
        .def_readwrite("_n_cam_use", &StereoMatch::_n_cam_use)
        .def_readwrite("_objID_match_list", &StereoMatch::_objID_match_list)
        .def_readwrite("_error_list", &StereoMatch::_error_list)
        .def("to_dict", [](StereoMatch const& self){
            CamList cam_list(self._cam_list);
            return py::dict(
                "_param"_a=self._param,
                "cam_list (no_access)"_a=cam_list, 
                "_n_cam_use"_a=self._n_cam_use, 
                "_objID_match_list"_a=self._objID_match_list, 
                "_error_list"_a=self._error_list
            );
        })
        .doc() = "StereoMatch class";
        
}

