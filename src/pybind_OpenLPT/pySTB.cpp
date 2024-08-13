
void init_STB(py::module &m) 
{
    py::class_<STB<Tracer3D>>(m, "STBTracer")
        .def(py::init<int, int, float, double, int, std::string const&, CamList const&, AxisLimit const&, std::string const&>(), py::arg("frame_start"), py::arg("frame_end"), py::arg("fps"), py::arg("vx_to_mm"), py::arg("n_thread"), py::arg("output_folder"), py::arg("cam_list"), py::arg("axis_limit"), py::arg("file"))
        .def("calibrateOTF", &STB<Tracer3D>::calibrateOTF)
        .def("processFrame", [](STB<Tracer3D>& self, int frame_id, std::vector<Image> const& img_list, bool is_update_img){
            std::vector<Image> img_list_copy(img_list);
            self.processFrame(frame_id, img_list_copy, is_update_img);  
            return img_list_copy;
        }, py::arg("frame_id"), py::arg("img_list"), py::arg("is_update_img")=false)
        .def(py::init<const STB<Tracer3D>&>())
        .def("loadTracks", &STB<Tracer3D>::loadTracks)
        .def("saveTracks", &STB<Tracer3D>::saveTracks)
        .def("saveTracksAll", &STB<Tracer3D>::saveTracksAll)
        .def("getObjParam", &STB<Tracer3D>::getObjParam)
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