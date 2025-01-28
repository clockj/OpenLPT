
void init_Track(py::module &m) 
{
    py::class_<TrackPredParam>(m, "TrackPredParam")
        .def(py::init<>())
        .def_readwrite("type", &TrackPredParam::type)
        .def_readwrite("param", &TrackPredParam::param)
        .def("to_dict", [](TrackPredParam const& self){
            return py::dict("type"_a=self.type, "param"_a=self.param);
        })
        .doc() = "TrackPredParam struct";

    py::class_<Track<Tracer3D>,  std::unique_ptr<Track<Tracer3D>>>(m, "TracerTrack")
        .def(py::init<>())
        .def(py::init<TrackPredParam const&>())
        .def(py::init<Tracer3D const&, int>())
        .def(py::init<Tracer3D const&, int, TrackPredParam const&>())
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
        .def_readwrite("_track_pred_param", &Track<Tracer3D>::_track_pred_param)
        .def("to_dict", [](Track<Tracer3D> const& self){
            return py::dict(
                "_obj3d_list"_a=self._obj3d_list, "_t_list"_a=self._t_list, "_n_obj3d"_a=self._n_obj3d, "_active"_a=self._active, "_track_pred_param"_a=self._track_pred_param
            );
        })
        .doc() = "TracerTrack class";
}
