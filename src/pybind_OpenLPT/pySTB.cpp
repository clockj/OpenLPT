#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <iostream>
#include <string>
#include <fstream>

#include "ObjectInfo.h"
#include "ObjectFinder.h"
#include "StereoMatch.h"
#include "OTF.h"
#include "Shake.h"
#include "IPR.h"
#include "PredField.h"
#include "Track.h"
#include "STB.h"

namespace py = pybind11;
using namespace pybind11::literals;

// Submodule: object
void init_ObjectInfo(py::module &m) 
{
    py::class_<Object2D>(m, "Object2D")
        .def(py::init<>())
        .def(py::init<Object2D const&>())
        .def(py::init<Pt2D const&>())
        .def_readwrite("_pt_center", &Object2D::_pt_center)
        .def("to_dict", [](Object2D const& self){
            return py::dict(
                "_pt_center"_a=self._pt_center
            );
        })
        .doc() = "Object2D class";

    py::class_<Tracer2D, Object2D>(m, "Tracer2D")
        .def(py::init<>())
        .def(py::init<Tracer2D const&>())
        .def(py::init<Pt2D const&>())
        .def_readwrite("_r_px", &Tracer2D::_r_px)
        .def("to_dict", [](Tracer2D const& self){
            return py::dict(
                "_pt_center"_a=self._pt_center, 
                "_r_px"_a=self._r_px
            );
        })
        .doc() = "Tracer2D class";

    py::class_<Object3D>(m, "Object3D")
        .def(py::init<>())
        .def(py::init<Object3D const&>())
        .def(py::init<Pt3D const&>())
        .def_readwrite("_pt_center", &Object3D::_pt_center)
        .def_readwrite("_is_tracked", &Object3D::_is_tracked)
        .def("to_dict", [](Object3D const& self){
            return py::dict(
                "_pt_center"_a=self._pt_center, 
                "_is_tracked"_a=self._is_tracked
            );
        })
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
        .def("to_dict", [](Tracer3D const& self){
            return py::dict(
                "_pt_center"_a=self._pt_center, 
                "_is_tracked"_a=self._is_tracked, 
                "_n_2d"_a=self._n_2d, 
                "_error"_a=self._error, 
                "_camid_list"_a=self._camid_list, 
                "_tr2d_list"_a=self._tr2d_list
            );
        })
        .doc() = "Tracer3D class";
}

void init_ObjectFinder(py::module& m)
{
    py::class_<ObjectFinder2D>(m, "ObjectFinder2D")
        .def(py::init<>())
        .def("findObject2D", [](ObjectFinder2D& self, Image const& img, std::vector<double> const& properties){
            std::vector<Tracer2D> obj2d_list;
            self.findObject2D(obj2d_list, img, properties);
            return obj2d_list;
        })
        .def("findObject2D", [](ObjectFinder2D& self, Image const& img, std::vector<double> const& properties, PixelRange const& region){
            std::vector<Tracer2D> obj2d_list;
            self.findObject2D(obj2d_list, img, properties, region);
            return obj2d_list;
        })
        .def("to_dict", [](ObjectFinder2D const& self){
            return py::dict();
        })
        .doc() = "ObjectFinder2D class";
}


// Submodule: stb
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

void init_OTF(py::module& m)
{
    py::class_<OTFParam>(m, "OTFParam")
        .def(py::init<>())
        .def_readwrite("dx", &OTFParam::dx)
        .def_readwrite("dy", &OTFParam::dy)
        .def_readwrite("dz", &OTFParam::dz)
        .def_readwrite("n_cam", &OTFParam::n_cam)
        .def_readwrite("nx", &OTFParam::nx)
        .def_readwrite("ny", &OTFParam::ny)
        .def_readwrite("nz", &OTFParam::nz)
        .def_readwrite("n_grid", &OTFParam::n_grid)
        .def_readwrite("a", &OTFParam::a)
        .def_readwrite("b", &OTFParam::b)
        .def_readwrite("c", &OTFParam::c)
        .def_readwrite("alpha", &OTFParam::alpha)
        .def_readwrite("boundary", &OTFParam::boundary)
        .def_readwrite("grid_x", &OTFParam::grid_x)
        .def_readwrite("grid_y", &OTFParam::grid_y)
        .def_readwrite("grid_z", &OTFParam::grid_z)
        .def("to_dict", [](OTFParam const& self){
            return py::dict(
                "dx"_a=self.dx, "dy"_a=self.dy, "dz"_a=self.dz, 
                "n_cam"_a=self.n_cam, "nx"_a=self.nx, "ny"_a=self.ny, "nz"_a=self.nz, "n_grid"_a=self.n_grid, 
                "a"_a=self.a, "b"_a=self.b, "c"_a=self.c, "alpha"_a=self.alpha, 
                "boundary"_a=self.boundary, "grid_x"_a=self.grid_x, "grid_y"_a=self.grid_y, "grid_z"_a=self.grid_z
            );
        })
        .doc() = "OTFParam struct";

    py::class_<OTF>(m, "OTF")
        .def(py::init<>())
        .def(py::init<const OTF&>())
        .def(py::init<int, int, int, int, AxisLimit const&>(), py::arg("n_cam"), py::arg("nx"), py::arg("ny"), py::arg("nz"), py::arg("boundary"))
        .def(py::init<std::string>())
        .def("loadParam", (void (OTF::*)(int, int, int, int, AxisLimit const&)) &OTF::loadParam, py::arg("n_cam"), py::arg("nx"), py::arg("ny"), py::arg("nz"), py::arg("boundary"))
        .def("loadParam", (void (OTF::*)(std::string)) &OTF::loadParam)
        .def("getOTFParam", &OTF::getOTFParam)
        .def_readwrite("_param", &OTF::_param)
        .def("to_dict", [](OTF const& self){
            return py::dict(
                "_param"_a=self._param
            );
        })
        .doc() = "OTF class";
}

void init_Shake(py::module& m)
{
    py::class_<Shake>(m, "Shake")
        .def(py::init<CamList const&, double, double, int, int>(), py::arg("cam_list"), py::arg("shake_width"), py::arg("score_min")=0.1, py::arg("n_loop")=4, py::arg("n_thread")=0)
        .def("runShake", [](Shake& self, std::vector<Tracer3D>const& obj3d_list, OTF const& otf, std::vector<Image> const& imgOrig_list, bool tri_only){
            std::vector<Tracer3D> tr3d_list_shake(obj3d_list);
            self.runShake(tr3d_list_shake, otf, imgOrig_list, tri_only);
            return tr3d_list_shake;
        }, py::arg("obj3d_list"), py::arg("otf"), py::arg("imgOrig_list"), py::arg("tri_only")=false)
        .def_readwrite("_imgRes_list", &Shake::_imgRes_list)
        .def_readwrite("_score_list", &Shake::_score_list)
        .def_readwrite("_objID_remove", &Shake::_objID_remove)
        .def_readwrite("_objID_keep", &Shake::_objID_keep)
        .def_readwrite("_n_ghost", &Shake::_n_ghost)
        .def("to_dict", [](Shake const& self){
            return py::dict(
                "_imgRes_list"_a=self._imgRes_list, "_score_list"_a=self._score_list, "_objID_remove"_a=self._objID_remove, "_objID_keep"_a=self._objID_keep, "_n_ghost"_a=self._n_ghost
            );
        })
        .doc() = "Shake class";
}

void init_IPR(py::module &m) 
{
    py::class_<IPRParam>(m, "IPRParam")
        .def(py::init<>())
        .def_readwrite("tri_only", &IPRParam::tri_only)
        .def_readwrite("n_thread", &IPRParam::n_thread)
        .def_readwrite("n_loop_ipr", &IPRParam::n_loop_ipr)
        .def_readwrite("n_loop_ipr_reduced", &IPRParam::n_loop_ipr_reduced)
        .def_readwrite("n_obj2d_max", &IPRParam::n_obj2d_max)
        .def_readwrite("tol_2d", &IPRParam::tol_2d)
        .def_readwrite("tol_3d", &IPRParam::tol_3d)
        .def_readwrite("check_id", &IPRParam::check_id)
        .def_readwrite("check_radius", &IPRParam::check_radius)
        .def_readwrite("n_loop_shake", &IPRParam::n_loop_shake)
        .def_readwrite("shake_width", &IPRParam::shake_width)
        .def_readwrite("ghost_threshold", &IPRParam::ghost_threshold)
        .def("to_dict", [](IPRParam const& self){
            return py::dict(
                "tri_only"_a=self.tri_only, "n_thread"_a=self.n_thread, "n_loop_ipr"_a=self.n_loop_ipr, "n_loop_ipr_reduced"_a=self.n_loop_ipr_reduced,
                "n_obj2d_max"_a=self.n_obj2d_max, "tol_2d"_a=self.tol_2d, "tol_3d"_a=self.tol_3d, "check_id"_a=self.check_id, "check_radius"_a=self.check_radius,
                "n_loop_shake"_a=self.n_loop_shake, "shake_width"_a=self.shake_width, "ghost_threshold"_a=self.ghost_threshold
            );
        })
        .doc() = "IPRParam struct";

    py::class_<IPR>(m, "IPR")
        .def(py::init<CamList&, std::vector<Image> const&, IPRParam const&>())
        .def("runIPR", [](IPR& self, std::vector<double> const& tr2d_properties, OTF const& otf, int n_reduced){
            std::vector<Tracer3D> tr3d_list_all;
            self.runIPR(tr3d_list_all, tr2d_properties, otf, n_reduced);
            return tr3d_list_all;
        
        }, py::arg("tr2d_properties"), py::arg("otf"), py::arg("n_reduced") = 0)
        .def("reducedCamLoop", [](IPR& self, std::vector<double> const& tr2d_properties, OTF const& otf, std::vector<int> const& cam_id, int n_cam){
            std::vector<Tracer3D> tr3d_list_all;
            self.reducedCamLoop(tr3d_list_all, tr2d_properties, otf, cam_id, n_cam);
            return tr3d_list_all;
        })
        .def("createCamID", [](IPR& self, std::vector<int> cam_id, int id, int n_rest){
            std::deque<std::vector<int>> cam_id_all;
            self.createCamID(cam_id_all, cam_id, id, n_rest);
            return cam_id_all;
        })
        .def("saveObjInfo", [](IPR& self, std::string const& filename, std::vector<Tracer3D> const& tr3d_list_all){
            self.saveObjInfo(filename, tr3d_list_all);
        })
        .def_readwrite("_param", &IPR::_param)
        .def_readwrite("_imgRes_list", &IPR::_imgRes_list)
        .def("to_dict", [](IPR const& self){
            return py::dict(
                "_param"_a=self._param, "_imgRes_list"_a=self._imgRes_list
            );
        })  
        .doc() = "IPR class";
}

void init_PredField(py::module &m) 
{
    py::class_<PFParam>(m, "PFParam")
        .def(py::init<>())
        .def_readwrite("limit", &PFParam::limit)
        .def_readwrite("nx", &PFParam::nx)
        .def_readwrite("ny", &PFParam::ny)
        .def_readwrite("nz", &PFParam::nz)
        .def_readwrite("r", &PFParam::r)
        .def_readwrite("nBin_x", &PFParam::nBin_x)
        .def_readwrite("nBin_y", &PFParam::nBin_y)
        .def_readwrite("nBin_z", &PFParam::nBin_z)
        .def("to_dict", [](PFParam const& self){
            return py::dict(
                "limit"_a=self.limit, "nx"_a=self.nx, "ny"_a=self.ny, "nz"_a=self.nz, "r"_a=self.r, "nBin_x"_a=self.nBin_x, "nBin_y"_a=self.nBin_y, "nBin_z"_a=self.nBin_z
            );
        })
        .doc() = "PFParam struct";

    py::class_<PredField>(m, "PredField")
        .def(py::init<PredField const&>())
        .def(py::init<PFParam const&, std::vector<Tracer3D> const&, std::vector<Tracer3D> const&>())
        .def(py::init<PFParam const&, std::vector<Pt3D> const&, std::vector<Pt3D> const&>())
        .def(py::init<PFParam const&, Matrix<double> const&>())
        .def("getDisp", &PredField::getDisp)
        .def("saveDispField", &PredField::saveDispField)
        .def_readwrite("_param", &PredField::_param)
        .def_readwrite("_disp_field", &PredField::_disp_field)
        .def("to_dict", [](PredField const& self){
            return py::dict(
                "_param"_a=self._param, "_disp_field"_a=self._disp_field
            );
        })
        .doc() = "PredField class";
}

void init_Track(py::module &m) 
{
    py::class_<Track<Tracer3D>,  std::unique_ptr<Track<Tracer3D>>>(m, "TracerTrack")
        .def(py::init<>())
        .def(py::init<Tracer3D const&, int>())
        .def(py::init<Track<Tracer3D> const&>())
        .def("addNext", (void (Track<Tracer3D>::*)(Tracer3D const&, int)) &Track<Tracer3D>::addNext)
        .def("addNext", (void (Track<Tracer3D>::*)(Track<Tracer3D> const&)) &Track<Tracer3D>::addNext)
        .def("predictNext", &Track<Tracer3D>::predictNext)
        .def("saveTrack", [](Track<Tracer3D>& self, std::string& file, int track_id, float fps, int n_cam_all){
            std::ofstream output(file, std::ios::app);
            self.saveTrack(output, track_id, fps, n_cam_all);
            output.close();
        }, py::arg("file"), py::arg("track_id"), py::arg("fps")=1, py::arg("n_cam_all")=0)
        .def_readwrite("_obj3d_list", &Track<Tracer3D>::_obj3d_list)
        .def_readwrite("_t_list", &Track<Tracer3D>::_t_list)
        .def_readwrite("_n_obj3d", &Track<Tracer3D>::_n_obj3d)
        .def_readwrite("_active", &Track<Tracer3D>::_active)
        .def("to_dict", [](Track<Tracer3D> const& self){
            return py::dict(
                "_obj3d_list"_a=self._obj3d_list, "_t_list"_a=self._t_list, "_n_obj3d"_a=self._n_obj3d, "_active"_a=self._active
            );
        })
        .doc() = "TracerTrack class";
}

void init_STB(py::module &m) 
{
    py::class_<STB<Tracer3D>>(m, "STBTracer")
        .def(py::init<int, int, float, double, int, std::string const&, CamList const&, AxisLimit const&, std::string const&>(), py::arg("frame_start"), py::arg("frame_end"), py::arg("fps"), py::arg("vx_to_mm"), py::arg("n_thread"), py::arg("output_folder"), py::arg("cam_list"), py::arg("axis_limit"), py::arg("file"))
        .def("calibrateOTF", &STB<Tracer3D>::calibrateOTF)
        .def("processFrame", &STB<Tracer3D>::processFrame, py::arg("frame_id"), py::arg("img_list"), py::arg("is_update_img")=false)
        .def("saveTracks", &STB<Tracer3D>::saveTracks)
        .def("saveTracksAll", &STB<Tracer3D>::saveTracksAll)
        .def_readwrite("_ipr_matched", &STB<Tracer3D>::_ipr_matched)
        .def_readwrite("_short_track_active", &STB<Tracer3D>::_short_track_active)
        .def_readwrite("_long_track_active", &STB<Tracer3D>::_long_track_active)
        .def_readwrite("_long_track_inactive", &STB<Tracer3D>::_long_track_inactive)
        .def_readwrite("_exit_track", &STB<Tracer3D>::_exit_track)
        .def("to_dict", [](STB<Tracer3D> const& self){
            return py::dict(
                "_ipr_matched"_a=self._ipr_matched, "_short_track_active"_a=self._short_track_active, "_long_track_active"_a=self._long_track_active, "_long_track_inactive"_a=self._long_track_inactive, "_exit_track"_a=self._exit_track
            );
        })
        .doc() = "STBTracer class";
}