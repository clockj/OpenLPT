
void init_Shake(py::module& m)
{
    py::class_<ImgAugList>(m, "ImgAugList")
        .def(py::init<>())
        .def_readwrite("img_list", &ImgAugList::img_list)
        .def_readwrite("region_list", &ImgAugList::region_list)
        .def("to_dict", [](ImgAugList const& self){
            return py::dict(
                "img_list"_a=self.img_list, "region_list"_a=self.region_list
            );
        })
        .doc() = "ImgAugList struct";

    py::class_<Shake>(m, "Shake")
        .def(py::init<CamList const&, double, double, double, int, int>(), py::arg("cam_list"), py::arg("shake_width"), py::arg("tol_3d"), py::arg("score_min")=0.1, py::arg("n_loop")=4, py::arg("n_thread")=0)
        .def("runShake", [](Shake& self, std::vector<Tracer3D>const& obj3d_list, OTF const& otf, std::vector<Image> const& imgOrig_list, bool tri_only){
            std::vector<Tracer3D> tr3d_list_shake(obj3d_list);
            self.runShake(tr3d_list_shake, otf, imgOrig_list, tri_only);
            return tr3d_list_shake;
        }, py::arg("obj3d_list"), py::arg("otf"), py::arg("imgOrig_list"), py::arg("tri_only")=false)
        .def_readwrite("_imgRes_list", &Shake::_imgRes_list)
        .def_readwrite("_score_list", &Shake::_score_list)
        .def_readwrite("_is_ghost", &Shake::_is_ghost)
        .def_readwrite("_is_repeated", &Shake::_is_repeated)
        .def_readwrite("_n_ghost", &Shake::_n_ghost)
        .def_readwrite("_n_repeated", &Shake::_n_repeated)
        .def("to_dict", [](Shake const& self){
            return py::dict(
                "_imgRes_list"_a=self._imgRes_list, "_score_list"_a=self._score_list, "_is_ghost"_a=self._is_ghost, "_is_repeated"_a=self._is_repeated, "_n_ghost"_a=self._n_ghost, "_n_repeated"_a=self._n_repeated
            );
        })
        .doc() = "Shake class";
}
