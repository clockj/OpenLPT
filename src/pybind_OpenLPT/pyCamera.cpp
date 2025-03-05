
void init_Camera(py::module &m) 
{
    py::class_<PinholeParam>(m, "PinholeParam")
        .def(py::init<>())
        .def_readwrite("n_row", &PinholeParam::n_row)
        .def_readwrite("n_col", &PinholeParam::n_col)
        .def_readwrite("cam_mtx", &PinholeParam::cam_mtx)
        .def_readwrite("is_distorted", &PinholeParam::is_distorted)
        .def_readwrite("n_dist_coeff", &PinholeParam::n_dist_coeff)
        .def_readwrite("dist_coeff", &PinholeParam::dist_coeff)
        .def_readwrite("r_mtx", &PinholeParam::r_mtx)
        .def_readwrite("t_vec", &PinholeParam::t_vec)
        .def_readwrite("r_mtx_inv", &PinholeParam::r_mtx_inv)
        .def_readwrite("t_vec_inv", &PinholeParam::t_vec_inv)
        .def("to_dict", [](PinholeParam const& self){
            return py::dict(
                "n_row"_a=self.n_row, 
                "n_col"_a=self.n_col, 
                "cam_mtx"_a=self.cam_mtx, 
                "is_distorted"_a=self.is_distorted, 
                "n_dist_coeff"_a=self.n_dist_coeff, 
                "dist_coeff"_a=self.dist_coeff, 
                "r_mtx"_a=self.r_mtx, 
                "t_vec"_a=self.t_vec, 
                "r_mtx_inv"_a=self.r_mtx_inv, 
                "t_vec_inv"_a=self.t_vec_inv
            );
        })
        .doc() = "PinholeParam struct";

    py::enum_<RefPlane>(m, "RefPlane")
        .value("REF_X", RefPlane::REF_X)
        .value("REF_Y", RefPlane::REF_Y)
        .value("REF_Z", RefPlane::REF_Z)
        .export_values();

    py::class_<PolyParam>(m, "PolyParam")
        .def(py::init<>())
        .def_readwrite("n_row", &PolyParam::n_row)
        .def_readwrite("n_col", &PolyParam::n_col)
        .def_readwrite("ref_plane", &PolyParam::ref_plane)
        .def_readwrite("plane", &PolyParam::plane)
        .def_readwrite("n_coeff", &PolyParam::n_coeff)
        .def_readwrite("u_coeffs", &PolyParam::u_coeffs)
        .def_readwrite("du_coeffs", &PolyParam::du_coeffs)
        .def_readwrite("v_coeffs", &PolyParam::v_coeffs)
        .def_readwrite("dv_coeffs", &PolyParam::dv_coeffs)
        .def("to_dict", [](PolyParam const& self){
            return py::dict(
                "n_row"_a=self.n_row, 
                "n_col"_a=self.n_col, 
                "ref_plane"_a=self.ref_plane, 
                "plane"_a=self.plane, 
                "n_coeff"_a=self.n_coeff, 
                "u_coeffs"_a=self.u_coeffs, 
                "du_coeffs"_a=self.du_coeffs, 
                "v_coeffs"_a=self.v_coeffs, 
                "dv_coeffs"_a=self.dv_coeffs
            );
        })
        .doc() = "PolyParam struct";
    
    py::class_<PinPlateParam, PinholeParam>(m, "PinPlateParam")
        .def(py::init<>())
        .def_readwrite("plane", &PinPlateParam::plane)
        .def_readwrite("pt3d_closest", &PinPlateParam::pt3d_closest)
        .def_readwrite("refract_array", &PinPlateParam::refract_array)
        .def_readwrite("w_array", &PinPlateParam::w_array)
        .def_readwrite("n_plate", &PinPlateParam::n_plate)
        .def_readwrite("proj_tol", &PinPlateParam::proj_tol)
        .def_readwrite("proj_nmax", &PinPlateParam::proj_nmax)
        .def_readwrite("lr", &PinPlateParam::lr)
        .def_readwrite("refract_ratio_max", &PinPlateParam::refract_ratio_max)
        .def("to_dict", [](PinPlateParam const& self){
            return py::dict(
                "n_row"_a=self.n_row, 
                "n_col"_a=self.n_col, 
                "cam_mtx"_a=self.cam_mtx, 
                "is_distorted"_a=self.is_distorted, 
                "n_dist_coeff"_a=self.n_dist_coeff, 
                "dist_coeff"_a=self.dist_coeff, 
                "r_mtx"_a=self.r_mtx, 
                "t_vec"_a=self.t_vec, 
                "r_mtx_inv"_a=self.r_mtx_inv, 
                "t_vec_inv"_a=self.t_vec_inv, 
                "plane"_a=self.plane, 
                "pt3d_closest"_a=self.pt3d_closest, 
                "refract_array"_a=self.refract_array, 
                "w_array"_a=self.w_array, 
                "n_plate"_a=self.n_plate, 
                "proj_tol"_a=self.proj_tol, 
                "proj_nmax"_a=self.proj_nmax,
                "lr"_a=self.lr,
                "refract_ratio_max"_a=self.refract_ratio_max
            );
        })
        .doc() = "PinPlateParam struct";

    py::enum_<CameraType>(m, "CameraType")
        .value("PINHOLE", CameraType::PINHOLE)
        .value("POLYNOMIAL", CameraType::POLYNOMIAL)
        .value("PINPLATE", CameraType::PINPLATE)
        .export_values();

    py::class_<Camera>(m, "Camera")
        .def_readwrite("_type", &Camera::_type)
        .def_readwrite("_pinhole_param", &Camera::_pinhole_param)
        .def_readwrite("_poly_param", &Camera::_poly_param)
        .def_readwrite("_pinplate_param", &Camera::_pinplate_param)
        .def(py::init<>())
        .def(py::init<const Camera&>())
        .def(py::init<std::string>())
        .def("loadParameters", [](Camera &self, std::string filename) {
            self.loadParameters(filename);
        })
        .def("updatePolyDuDv", &Camera::updatePolyDuDv)
        .def("updatePt3dClosest", &Camera::updatePt3dClosest)
        .def("saveParameters", &Camera::saveParameters)
        .def("rmtxTorvec", &Camera::rmtxTorvec)
        .def("getNRow", &Camera::getNRow)
        .def("getNCol", &Camera::getNCol)
        .def("project", [](Camera const& self, Pt3D const& pt3d, bool is_print_detail){
            return self.project(pt3d, is_print_detail);
        }, py::arg("pt3d"), py::arg("is_print_detail")=false)
        .def("project", [](Camera const& self, std::vector<Pt3D> const& pt3d_list){
            std::vector<Pt2D> pt2d_list(pt3d_list.size());
            #pragma omp parallel for
            for (int i = 0; i < pt3d_list.size(); i++)
            {
                pt2d_list[i] = self.project(pt3d_list[i]);
            }
            return pt2d_list;
        })
        .def("worldToUndistImg", &Camera::worldToUndistImg)
        .def("worldToUndistImg", [](Camera const& self, std::vector<Pt3D> const& pt3d_list, PinholeParam const& param){
            int npts = pt3d_list.size();
            std::vector<Pt2D> pt2d_list(npts);
            #pragma omp parallel for
            for (int i = 0; i < npts; i++)
            {
                pt2d_list[i] = self.worldToUndistImg(pt3d_list[i], param);
            }
            return pt2d_list;
        })
        .def("refractPlate", &Camera::refractPlate)
        .def("distort", &Camera::distort)
        .def("distort", [](Camera const& self, std::vector<Pt2D> const& pt2d_list, PinholeParam const& param){
            int npts = pt2d_list.size();
            std::vector<Pt2D> pt2d_list_dist(npts);
            #pragma omp parallel for
            for (int i = 0; i < npts; i++)
            {
                pt2d_list_dist[i] = self.distort(pt2d_list[i], param);
            }
            return pt2d_list_dist;
        })
        .def("polyProject", &Camera::polyProject)
        .def("lineOfSight", &Camera::lineOfSight)
        .def("lineOfSight", [](Camera const& self, std::vector<Pt2D> const& pt2d_list){
            std::vector<Line3D> line_list(pt2d_list.size());
            #pragma omp parallel for
            for (int i = 0; i < pt2d_list.size(); i++)
            {
                line_list[i] = self.lineOfSight(pt2d_list[i]);
            }
            return line_list;
        })
        .def("undistort", &Camera::undistort)
        .def("undistort", [](Camera const& self, std::vector<Pt2D> const& pt2d_list, PinholeParam const& param){
            int npts = pt2d_list.size();
            std::vector<Pt2D> pt2d_list_undist(npts);
            #pragma omp parallel for
            for (int i = 0; i < npts; i++)
            {
                pt2d_list_undist[i] = self.undistort(pt2d_list[i], param);
            }
            return pt2d_list_undist;
        })
        .def("pinholeLine", &Camera::pinholeLine)
        .def("polyImgToWorld", &Camera::polyImgToWorld)
        .def("polyLineOfSight", &Camera::polyLineOfSight)
        .def("pinplateLine", &Camera::pinplateLine)
        .def("to_dict", [](Camera const& self){
            return py::dict(
                "_type"_a=self._type, 
                "_pinhole_param"_a=self._pinhole_param, 
                "_poly_param"_a=self._poly_param,
                "_pinplate_param"_a=self._pinplate_param
            );
        })
        .doc() = "Camera class";

    py::class_<CamList>(m, "CamList")
        .def(py::init<>())
        .def(py::init([](std::vector<Camera> const& cam_list, std::vector<int> const& intensity_max, std::vector<int> const& useid_list){
            CamList cam_list_all;
            cam_list_all.cam_list = cam_list;
            cam_list_all.intensity_max = intensity_max;
            cam_list_all.useid_list = useid_list;
            return cam_list_all;
        }))
        .def_readwrite("cam_list", &CamList::cam_list)
        .def_readwrite("intensity_max", &CamList::intensity_max)
        .def_readwrite("useid_list", &CamList::useid_list)
        .def("to_dict", [](CamList const& self){
            return py::dict(
                "cam_list"_a=self.cam_list, 
                "intensity_max"_a=self.intensity_max, 
                "useid_list"_a=self.useid_list
            );
        })
        .doc() = "CamList struct";
}