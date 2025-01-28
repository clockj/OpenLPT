void init_KalmanFilter(py::module &m) 
{
    py::class_<KalmanFilter>(m, "KalmanFilter")
        .def(py::init<>())
        .def(py::init<Matrix<double> const&, Matrix<double> const&, Matrix<double> const&, Matrix<double> const&, Matrix<double> const&, Matrix<double> const&>())
        .def("predict", &KalmanFilter::predict)
        .def("update", &KalmanFilter::update)
        .def_readwrite("_F", &KalmanFilter::_F)
        .def_readwrite("_H", &KalmanFilter::_H)
        .def_readwrite("_Q", &KalmanFilter::_Q)
        .def_readwrite("_R", &KalmanFilter::_R)
        .def_readwrite("_x", &KalmanFilter::_x)
        .def_readwrite("_P", &KalmanFilter::_P)
        .def("to_dict", [](KalmanFilter const& self){
            return py::dict(
                "_F"_a=self._F, 
                "_H"_a=self._H, 
                "_Q"_a=self._Q, 
                "_R"_a=self._R, 
                "_x"_a=self._x, 
                "_P"_a=self._P
            );
        })
        .doc() = "KalmanFilter class";
}