
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
        .def("findTracer2D", [](ObjectFinder2D& self, Image const& img, std::vector<double> const& properties){
            std::vector<Tracer2D> obj2d_list;
            self.findObject2D(obj2d_list, img, properties);
            return obj2d_list;
        })
        .def("findTracer2D", [](ObjectFinder2D& self, Image const& img, std::vector<double> const& properties, PixelRange const& region){
            std::vector<Tracer2D> obj2d_list;
            self.findObject2D(obj2d_list, img, properties, region);
            return obj2d_list;
        })
        .def("to_dict", [](ObjectFinder2D const& self){
            return py::dict();
        })
        .doc() = "ObjectFinder2D class";
}
