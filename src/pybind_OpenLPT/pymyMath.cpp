
void init_myMath(py::module& m)
{
    // Sort index
    m.def("sortID", [](std::vector<double> const& nums){
        std::vector<int> idx(nums.size());
        myMATH::sortID(idx, nums);
        return idx;
    }, "Sort ID for vector<double>");
    
    // Calculate median
    m.def("getMedian", &myMATH::getMedian<double>, "Get median for vector<double>");
    
    m.def("isOutlier", [](std::vector<double> const& nums){
        std::vector<int> judge(nums.size());
        myMATH::isOutlier(judge, nums);
        return judge;
    }, "Check if value is an outlier for vector<double>");
    
    m.def("linspace", &myMATH::linspace, "Create a linearly spaced vector<double>");
    
    m.def("triLinearInterp", &myMATH::triLinearInterp, "Perform trilinear interpolation for vector<double>");
    
    m.def("createUnitVector", [](Pt3D const& pt1, Pt3D const& pt2){
        return myMATH::createUnitVector(pt1, pt2);
    }, "Create a unit vector for two points");
    m.def("createUnitVector", [](Pt2D const& pt1, Pt2D const& pt2){
        return myMATH::createUnitVector(pt1, pt2);
    }, "Create a unit vector for three points");
    
    m.def("dot", [](Pt3D const& pt1, Pt3D const& pt2){
        return myMATH::dot(pt1, pt2);
    }, "Dot product for two 3D points");
    m.def("dot", [](Pt2D const& pt1, Pt2D const& pt2){
        return myMATH::dot(pt1, pt2);
    }, "Dot product for two 2D points");

    m.def("dist2", [](Pt3D const& pt1, Pt3D const& pt2){
        return myMATH::dist2(pt1, pt2);
    }, "Squared distance between two 3D points");
    m.def("dist2", [](Pt2D const& pt1, Pt2D const& pt2){
        return myMATH::dist2(pt1, pt2);
    }, "Squared distance between two 2D points");

    m.def("dist2", [](Pt3D const& pt, Line3D const& line){
        return myMATH::dist2(pt, line);
    }, "Squared distance between a 3D point and a 3D line");
    m.def("dist2", [](Pt2D const& pt, Line2D const& line){
        return myMATH::dist2(pt, line);
    }, "Squared distance between a 2D point and a 2D line");

    m.def("dist", [](Pt3D const& pt1, Pt3D const& pt2){
        return myMATH::dist(pt1, pt2);
    }, "Distance between two 3D points");
    m.def("dist", [](Pt2D const& pt1, Pt2D const& pt2){
        return myMATH::dist(pt1, pt2);
    }, "Distance between two 2D points");

    m.def("dist", [](Pt3D const& pt, Line3D const& line){
        return myMATH::dist(pt, line);
    }, "Distance between a 3D point and a 3D line");
    m.def("dist", [](Pt2D const& pt, Line2D const& line){
        return myMATH::dist(pt, line);
    }, "Distance between a 2D point and a 2D line");

    m.def("triangulation", [](std::vector<Line3D> const& line_of_sight_list){
        Pt3D pt_world(0, 0, 0);
        double error = 0;
        myMATH::triangulation(pt_world, error, line_of_sight_list);
        return std::make_pair(pt_world, error);
    }, "Triangulate 3D points from a list of 3D lines of sight");

    m.def("crossPoint", &myMATH::crossPoint, "Find the cross point of two 2D lines");

    m.def("eye", &myMATH::eye<double>, "Create an identity Matrix<double>");

    m.def("piecewiseProduct", &myMATH::piecewiseProduct<double>, "Piecewise product for two Matrix<double>");

    m.def("inverse", &myMATH::inverse<double>, "Inverse for Matrix<double>");

    m.def("trace", &myMATH::trace<double>, "Trace for Matrix<double>");

    m.def("isLocalMax", &myMATH::isLocalMax<double>, "Check if value is a local maximum for Matrix<double>");

    m.def("polyfit", [](std::vector<double> const& x, std::vector<double> const& y, int order){
        std::vector<double> coeff(order + 1);
        myMATH::polyfit(coeff, x, y, order);
        return coeff;
    }, "Polynomial fitting for vector<double>");
}