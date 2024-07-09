
void init_ImageIO(py::module &m) 
{
    py::class_<ImageParam>(m, "ImageParam")
        .def(py::init<>())
        .def_readwrite("n_row", &ImageParam::n_row)
        .def_readwrite("n_col", &ImageParam::n_col)
        .def_readwrite("bits_per_sample", &ImageParam::bits_per_sample)
        .def_readwrite("n_channel", &ImageParam::n_channel)
        .def("to_dict", [](ImageParam const& self){
            return py::dict(
                "n_row"_a=self.n_row, 
                "n_col"_a=self.n_col, 
                "bits_per_sample"_a=self.bits_per_sample, 
                "n_channel"_a=self.n_channel
            );
        })
        .doc() = "ImageParam struct";

    py::class_<ImageIO>(m, "ImageIO")
        .def(py::init<>())
        .def(py::init<ImageIO const&>())
        .def(py::init<std::string, std::string>())
        .def("loadImgPath", &ImageIO::loadImgPath)
        .def("loadImg", &ImageIO::loadImg)
        .def("saveImg", &ImageIO::saveImg)
        .def("setImgParam", &ImageIO::setImgParam)
        .def("getImgParam", &ImageIO::getImgParam)
        .def("getCurrImgID", &ImageIO::getCurrImgID)
        .def("to_dict", [](ImageIO const& self){
            return py::dict(
                "img_path (no_access)"_a=self.getImgPath(), 
                "img_param (no_access)"_a=self.getImgParam(),
                "curr_img_id (no_access)"_a=self.getCurrImgID()
            );
        })
        .doc() = "ImageIO class";
}