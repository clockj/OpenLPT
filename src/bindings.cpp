#include <pybind11/pybind11.h>
#include <iostream>
#include <streambuf>
#include <string>

namespace py = pybind11;

// STBCommons struct and constants
void init_STBCommons(py::module &);

// Submodule: math
void init_Matrix(py::module &);
void init_ImageIO(py::module &);
void init_Camera(py::module &);
void init_myMath(py::module &);


// Redirect std::cout to Python's sys.stdout
class PythonStreamRedirector : public std::streambuf 
{
public:
    PythonStreamRedirector() : default_buffer(std::cout.rdbuf()) 
    {
        std::cout.rdbuf(this);
    }

    ~PythonStreamRedirector() 
    {
        std::cout.rdbuf(default_buffer);
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
            std::string s(1, static_cast<char>(c));
            
            py::module_::import("sys").attr("stdout").attr("write")(s);

            return c;
        }
    }

    virtual std::streamsize xsputn(const char* s, std::streamsize n) override 
    {
        std::string str(s, n);
        if (!str.empty()) 
        {  
            py::module_::import("sys").attr("stdout").attr("write")(str);
        }
        return n;
    }

private:
    std::streambuf* default_buffer;

    void print_buffer() 
    {
        std::string buffer = std::string(pbase(), pptr() - pbase());
        py::module_::import("sys").attr("stdout").attr("write")(buffer);
        py::print(buffer);
        setp(pbase(), epptr());
    }
};


// Define the module
PYBIND11_MODULE(pyOpenLPT, m) 
{
    // Redirect std::cout to Python's sys.stdout
    py::class_<PythonStreamRedirector>(m, "PythonStreamRedirector")
        .def(py::init<>());

    
    // STBCommons struct and constants
    init_STBCommons(m);


    // Submodule: math
    py::module m_math = m.def_submodule("math", "Math module");
    init_Matrix(m_math);
    init_ImageIO(m_math);
    init_myMath(m_math);
    init_Camera(m_math);
    
}