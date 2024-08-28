
template <class T>
py::array_t<T> matrix_to_numpy(Matrix<T> const& mat) 
{
    // Define the shape of the array (rows, cols)
    size_t n_row = mat.getDimRow();
    size_t n_col = mat.getDimCol();
    std::vector<size_t> shape = {n_row, n_col};

    // Create a NumPy array with the same size and data type as the matrix
    py::array_t<T> array(shape);

    // Copy data from matrix to NumPy array
    auto buf = array.request();
    T* ptr = static_cast<T*>(buf.ptr);
    std::copy(mat.data(), mat.data() + (n_row*n_col), ptr);

    return array;
}


template <typename T>
Matrix<T> numpy_to_matrix(const py::array_t<T>& array) 
{
    auto buf = array.request();
    if (buf.ndim != 2) 
    {
        throw std::runtime_error("NumPy array must be 2-dimensional");
    }
    int rows = buf.shape[0];
    int cols = buf.shape[1];
    Matrix<T> mat(rows, cols, 0);
    T* ptr = static_cast<T*>(buf.ptr);
    mat.setData(ptr, rows * cols);
    return mat;
}


py::array_t<double> pts_to_numpy(std::vector<Pt3D> const& pt3d_list)
{
    size_t n_row = pt3d_list.size();
    size_t n_col = 3;
    std::vector<size_t> shape = {n_row, n_col};
    py::array_t<double> array(shape);
    auto buf = array.request();
    double* ptr = static_cast<double*>(buf.ptr);
    
    #pragma omp parallel for
    for (int i = 0; i < n_row; i++)
    {
        ptr[i*n_col] = pt3d_list[i][0];
        ptr[i*n_col+1] = pt3d_list[i][1];
        ptr[i*n_col+2] = pt3d_list[i][2];
    }

    return array;
}

py::array_t<double> pts_to_numpy(std::vector<Pt2D> const& pt2d_list)
{
    size_t n_row = pt2d_list.size();
    size_t n_col = 2;
    std::vector<size_t> shape = {n_row, n_col};
    py::array_t<double> array(shape);
    auto buf = array.request();
    double* ptr = static_cast<double*>(buf.ptr);

    #pragma omp parallel for
    for (int i = 0; i < n_row; i++)
    {
        ptr[i*n_col] = pt2d_list[i][0];
        ptr[i*n_col+1] = pt2d_list[i][1];
    }

    return array;
}


std::vector<Pt3D> numpy_to_pt3d(py::array_t<double> const& array)
{
    auto buf = array.request();
    if (buf.ndim != 2) 
    {
        throw std::runtime_error("NumPy array must be 2-dimensional");
    }
    int rows = buf.shape[0];
    int cols = buf.shape[1];
    if (cols != 3) 
    {
        throw std::runtime_error("NumPy array must have 3 columns");
    }
    double* ptr = static_cast<double*>(buf.ptr);
    std::vector<Pt3D> pt3d_list(rows);
    #pragma omp parallel for
    for (int i = 0; i < rows; i++)
    {
        pt3d_list[i] = Pt3D(ptr[i*cols], ptr[i*cols+1], ptr[i*cols+2]);
    }
    return pt3d_list;
}

std::vector<Pt2D> numpy_to_pt2d(py::array_t<double> const& array)
{
    auto buf = array.request();
    if (buf.ndim != 2) 
    {
        throw std::runtime_error("NumPy array must be 2-dimensional");
    }
    int rows = buf.shape[0];
    int cols = buf.shape[1];
    if (cols != 2) 
    {
        throw std::runtime_error("NumPy array must have 2 columns");
    }
    double* ptr = static_cast<double*>(buf.ptr);
    std::vector<Pt2D> pt2d_list(rows);
    #pragma omp parallel for
    for (int i = 0; i < rows; i++)
    {
        pt2d_list[i] = Pt2D(ptr[i*cols], ptr[i*cols+1]);
    }
    return pt2d_list;
}


void init_Matrix(py::module &m) 
{
    py::class_<Matrix<double>, std::unique_ptr<Matrix<double>>>(m, "Matrix")
        .def(py::init<>())
        .def(py::init<Matrix<double> const&>())
        .def(py::init<int, int, double>())
        .def(py::init<std::string>())
        .def("__setitem__", [](Matrix<double> &self, int id, float val) {self[id] = val;})
        .def("__setitem__", [](Matrix<double> &self, std::pair<int,int> index, float val) {
            self(index.first,index.second) = val;
        })
        .def("__getitem__", [](Matrix<double> const&self, int id) {return self[id];})
        .def("__getitem__", [](Matrix<double> const&self, std::pair<int,int> index) {
            return self(index.first,index.second);
        })
        .def("getRow", &Matrix<double>::getRow)
        .def("getCol", &Matrix<double>::getCol)
        .def("getDimRow", &Matrix<double>::getDimRow)
        .def("getDimCol", &Matrix<double>::getDimCol)
        .def("print", &Matrix<double>::print, py::arg("precision") = 3)
        .def("write", [](Matrix<double> &self, std::string file_name) {
            self.write(file_name);
        })
        .def("norm", &Matrix<double>::norm)
        .def("__eq__", &Matrix<double>::operator==)
        .def("__ne__", &Matrix<double>::operator!=)
        .def("__add__", [](Matrix<double> &self, Matrix<double> const& mtx) {
            return self + mtx;
        })
        .def("__iadd__", [](Matrix<double> &self, Matrix<double> const& mtx) {
            self += mtx;
            return self;
        })
        .def("__add__", [](Matrix<double> &self, double& delta) {
            return self + delta;
        })
        .def("__iadd__", [](Matrix<double> &self, double& delta) {
            self += delta;
            return self;
        })
        .def("__sub__", [](Matrix<double> &self, Matrix<double> const& mtx) {
            return self - mtx;
        })
        .def("__isub__", [](Matrix<double> &self, Matrix<double> const& mtx) {
            self -= mtx;
            return self;
        })
        .def("__sub__", [](Matrix<double> &self, double& delta) {
            return self - delta;
        })
        .def("__isub__", [](Matrix<double> &self, double& delta) {
            self -= delta;
            return self;
        })
        .def("__mul__", [](Matrix<double> &self, Matrix<double> const& mtx) {
            return self * mtx;
        })
        .def("__imul__", [](Matrix<double> &self, Matrix<double> const& mtx) {
            self *= mtx;
            return self;
        })
        .def("__mul__", [](Matrix<double> &self, double& ratio) {
            return self * ratio;
        })
        .def("__imul__", [](Matrix<double> &self, double& ratio) {
            self *= ratio;
            return self;
        })
        .def("__truediv__", &Matrix<double>::operator/)
        .def("__itruediv__", &Matrix<double>::operator/=)
        .def("transpose", &Matrix<double>::transpose)
        .def_property_readonly("T", &Matrix<double>::transpose) 
        .def("to_dict", [](Matrix<double> const& self){
            return py::dict(
                "data (no_access)"_a=matrix_to_numpy<double>(self)
            );
        })
        .doc() = "Matrix<double> class";
    
    py::class_<Pt3D, Matrix<double>>(m, "Pt3D")
        .def(py::init<>())
        .def(py::init<double, double, double>())
        .def(py::init<Pt3D const&>())
        .def(py::init<Matrix<double> const&>())
        .def(py::init<std::string>())
        .doc() = "Pt3D class";

    py::class_<Pt2D, Matrix<double>>(m, "Pt2D")
        .def(py::init<>())
        .def(py::init<double, double>())
        .def(py::init<Pt2D const&>())
        .def(py::init<Matrix<double> const&>())
        .def(py::init<std::string>())
        .doc() = "Pt2D class";

    py::class_<Line3D>(m, "Line3D")
        .def(py::init<>())
        .def(py::init<const Pt3D&, const Pt3D&>())
        .def_readwrite("pt", &Line3D::pt)
        .def_readwrite("unit_vector", &Line3D::unit_vector)
        .def("to_dict", [](Line3D const& self){
            return py::dict(
                "pt"_a=self.pt, 
                "unit_vector"_a=self.unit_vector
            );
        })
        .doc() = "Line3D struct";

    py::class_<Line2D>(m, "Line2D")
        .def(py::init<>())
        .def(py::init<const Pt2D&, const Pt2D&>())
        .def_readwrite("pt", &Line2D::pt)
        .def_readwrite("unit_vector", &Line2D::unit_vector)
        .def("to_dict", [](Line2D const& self){
            return py::dict(
                "pt"_a=self.pt, 
                "unit_vector"_a=self.unit_vector
            );
        })
        .doc() = "Line2D struct";

    py::class_<Image, Matrix<double>>(m, "Image")
        .def(py::init<>())
        .def(py::init<int, int, double>())
        .def(py::init<Image const&>())
        .def(py::init<Matrix<double> const&>())
        .def(py::init<std::string>())
        .doc() = "Image class";

    m.def("matrix_to_numpy", &matrix_to_numpy<double>, "Convert a Matrix<double> to a NumPy array");
    m.def("numpy_to_matrix", &numpy_to_matrix<double>, "Convert a NumPy array to a Matrix<double>");
    m.def("pts_to_numpy", [](std::vector<Pt2D> const& pt2d_list){
        return pts_to_numpy(pt2d_list);
    }, "Convert a list of Pt2D to a NumPy array");
    m.def("pts_to_numpy", [](std::vector<Pt3D> const& pt3d_list){
        return pts_to_numpy(pt3d_list);
    }, "Convert a list of Pt3D to a NumPy array");
    m.def("numpy_to_pt2d", &numpy_to_pt2d, "Convert a NumPy array to a list of Pt2D");
    m.def("numpy_to_pt3d", &numpy_to_pt3d, "Convert a NumPy array to a list of Pt3D");
}

