# Compile static lib to be used in the pybind11 module
# Math module
add_library(bindMatrix INTERFACE ${CMAKE_SOURCE_DIR}/src/srcMath/Matrix.hpp)
set_property(TARGET bindMatrix PROPERTY LINKER_LANGUAGE CXX)

add_library(bindImageIO STATIC ${CMAKE_SOURCE_DIR}/src/srcMath/ImageIO.cpp)
add_subdirectory("${CMAKE_HOME_DIRECTORY}/inc/libtiff")
target_link_libraries(bindImageIO PUBLIC tiff)
# target_link_libraries(bindImageIO PUBLIC bindMatrix tiff)

add_library(bindmyMath STATIC ${CMAKE_SOURCE_DIR}/src/srcMath/myMath.cpp)
# target_link_libraries(bindmyMath PUBLIC bindMatrix)

add_library(bindCamera STATIC ${CMAKE_SOURCE_DIR}/src/srcMath/Camera.cpp)
# target_link_libraries(bindCamera PUBLIC bindMatrix bindmyMath)


# Object module
add_library(bindObjectInfo STATIC ${CMAKE_SOURCE_DIR}/src/srcObject/ObjectInfo.cpp)
# target_link_libraries(bindObjectInfo PUBLIC bindMatrix)

add_library(bindObjectFinder INTERFACE ${CMAKE_SOURCE_DIR}/src/srcObject/ObjectFinder.hpp)
set_property(TARGET bindObjectFinder PROPERTY LINKER_LANGUAGE CXX)
# target_link_libraries(bindObjectFinder PUBLIC bindMatrix bindObjectInfo bindmyMath)


# STB module 
# Find openmp package
add_library(bindStereoMatch INTERFACE ${CMAKE_SOURCE_DIR}/src/srcSTB/StereoMatch.hpp)
set_property(TARGET bindStereoMatch PROPERTY LINKER_LANGUAGE CXX)

add_library(bindOTF STATIC ${CMAKE_SOURCE_DIR}/src/srcSTB/OTF.cpp)

add_library(bindShake STATIC ${CMAKE_SOURCE_DIR}/src/srcSTB/Shake.cpp)

add_library(bindIPR INTERFACE ${CMAKE_SOURCE_DIR}/src/srcSTB/IPR.hpp)
set_property(TARGET bindIPR PROPERTY LINKER_LANGUAGE CXX)

add_library(bindPredField INTERFACE ${CMAKE_SOURCE_DIR}/src/srcSTB/PredField.hpp)
set_property(TARGET bindPredField PROPERTY LINKER_LANGUAGE CXX)

add_library(bindTrack INTERFACE ${CMAKE_SOURCE_DIR}/src/srcSTB/Track.hpp)
set_property(TARGET bindTrack PROPERTY LINKER_LANGUAGE CXX)

add_library(bindSTB INTERFACE ${CMAKE_SOURCE_DIR}/src/srcSTB/STB.hpp)
set_property(TARGET bindSTB PROPERTY LINKER_LANGUAGE CXX)


# Create pybind11 module
set(BINDINGS_SRC
    ${CMAKE_SOURCE_DIR}/src/bindings.cpp
    ${CMAKE_SOURCE_DIR}/src/pybind_OpenLPT/pySTBCommons.cpp
    ${CMAKE_SOURCE_DIR}/src/pybind_OpenLPT/pyMatrix.cpp
    ${CMAKE_SOURCE_DIR}/src/pybind_OpenLPT/pyImageIO.cpp
    ${CMAKE_SOURCE_DIR}/src/pybind_OpenLPT/pymyMath.cpp
    ${CMAKE_SOURCE_DIR}/src/pybind_OpenLPT/pyCamera.cpp
    ${CMAKE_SOURCE_DIR}/src/pybind_OpenLPT/pySTB.cpp
)
set(BINDINGS_LIB
    bindMatrix
    bindImageIO
    bindmyMath
    bindCamera
    bindObjectInfo
    bindObjectFinder
    bindStereoMatch
    bindOTF
    bindShake
    bindIPR
    bindPredField
    bindTrack
    bindSTB
)
pybind11_add_module(pyOpenLPT ${BINDINGS_SRC})
target_link_libraries(pyOpenLPT PRIVATE ${BINDINGS_LIB})

