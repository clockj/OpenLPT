#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <iostream>
#include <string>
#include <tuple>
#include <sstream>
#include <fstream>

#include "STBCommons.h"
#include "ImageIO.h"
#include "Matrix.h"
#include "Camera.h"
#include "myMATH.h"
#include "ObjectInfo.h"
#include "ObjectFinder.h"
#include "StereoMatch.h"
#include "OTF.h"
#include "Shake.h"
#include "IPR.h"
#include "PredField.h"
#include "Track.h"
#include "STB.h"
#include "main.cpp"

namespace py = pybind11;
using namespace pybind11::literals;

#include "pySTBCommons.cpp"
#include "pyOpenLPT.cpp"

#include "pyMatrix.cpp"
#include "pyImageIO.cpp"
#include "pyCamera.cpp"
#include "pymyMath.cpp"

#include "pyObjectInfo.cpp"
#include "pyObjectFinder.cpp"

#include "pyStereoMatch.cpp"
#include "pyOTF.cpp"
#include "pyShake.cpp"
#include "pyIPR.cpp"
#include "pyPredField.cpp"
#include "pyTrack.cpp"
#include "pySTB.cpp"
