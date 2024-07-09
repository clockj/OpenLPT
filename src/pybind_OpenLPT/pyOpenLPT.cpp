
void init_OpenLPT(py::module &m) 
{
    m.def("run", &run, "run OpenLPT");
}
