#include <pybind11/pybind11.h>
#include <iostream>
#include <streambuf>
#include <string>

namespace py = pybind11;

// STBCommons struct and constants
void init_STBCommons(py::module &);
void init_OpenLPT(py::module &);

// Submodule: math
void init_Matrix(py::module &);
void init_ImageIO(py::module &);
void init_Camera(py::module &);
void init_myMath(py::module &);

// Submodule: object 
void init_ObjectInfo(py::module &);
void init_ObjectFinder(py::module &);

// Submodule: stb
void init_StereoMatch(py::module &);
void init_OTF(py::module &);
void init_Shake(py::module &);
void init_IPR(py::module &);
void init_PredField(py::module &);
void init_Track(py::module &);
void init_STB(py::module &);


// Redirect std::cout to Python's sys.stdout
class PythonStreamRedirector : public std::streambuf {
public:
    PythonStreamRedirector() : default_cout_buffer(std::cout.rdbuf()), default_cerr_buffer(std::cerr.rdbuf())
    {
        std::cout.rdbuf(this);
        std::cerr.rdbuf(this);
    }

    ~PythonStreamRedirector() 
    {
        std::cout.rdbuf(default_cout_buffer);
        std::cerr.rdbuf(default_cerr_buffer);
    }

protected:
    virtual int overflow(int c) override 
    {
        if (c == EOF) 
        {
            return !EOF;
        } 
        else 
        {
            char z = c;
            py::module_::import("sys").attr("stderr").attr("write")(std::string(1, z));
            return c;
        }
    }

    virtual std::streamsize xsputn(const char* s, std::streamsize n) override 
    {
        std::string str(s, n);
        py::module_::import("sys").attr("stderr").attr("write")(str);
        return n;
    }

private:
    std::streambuf* default_cout_buffer;
    std::streambuf* default_cerr_buffer;
}


// Define the module
PYBIND11_MODULE(pyOpenLPT, m) 
{
    // Redirect std::cout to Python's sys.stdout
    py::class_<PythonStreamRedirector>(m, "PythonStreamRedirector")
        .def(py::init<>());

    
    // STBCommons struct and constants
    init_STBCommons(m);
    init_OpenLPT(m);

    // Submodule: math
    py::module m_math = m.def_submodule("math", "Math module");
    init_Matrix(m_math);
    init_ImageIO(m_math);
    init_myMath(m_math);
    init_Camera(m_math);
    

    // Submodule: object
    py::module m_object = m.def_submodule("object", "Object module");
    init_ObjectInfo(m_object);
    init_ObjectFinder(m_object);


    // Submodule: stb
    py::module m_stb = m.def_submodule("stb", "STB module");
    init_StereoMatch(m_stb);
    init_OTF(m_stb);
    init_Shake(m_stb);
    init_IPR(m_stb);
    init_PredField(m_stb);
    init_Track(m_stb);
    init_STB(m_stb);
}