
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
        .def("saveParam", &OTF::saveParam)
        .def("getOTFParam", &OTF::getOTFParam)
        .def_readwrite("_param", &OTF::_param)
        .def("to_dict", [](OTF const& self){
            return py::dict(
                "_param"_a=self._param
            );
        })
        .doc() = "OTF class";
}
