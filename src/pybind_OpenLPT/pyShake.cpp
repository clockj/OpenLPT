
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
