
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
