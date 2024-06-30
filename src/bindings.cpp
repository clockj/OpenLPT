#include <pybind11/pybind11.h>
#include <iostream>
#include <streambuf>
#include <string>

namespace py = pybind11;

void init_Matrix(py::module &);

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


PYBIND11_MODULE(pyOpenLPT, m) 
{
    py::class_<PythonStreamRedirector>(m, "PythonStreamRedirector")
        .def(py::init<>());

    py::module m_math = m.def_submodule("math", "Math module");

    init_Matrix(m_math);
}