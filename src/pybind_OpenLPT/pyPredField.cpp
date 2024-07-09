
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
