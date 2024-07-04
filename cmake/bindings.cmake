# Compile static lib to be used in the pybind11 module
add_library(bindMatrix INTERFACE ${CMAKE_SOURCE_DIR}/src/srcMath/Matrix.hpp)
set_property(TARGET bindMatrix PROPERTY LINKER_LANGUAGE CXX)

add_library(bindImageIO STATIC ${CMAKE_SOURCE_DIR}/src/srcMath/ImageIO.cpp)
add_subdirectory("${CMAKE_HOME_DIRECTORY}/inc/libtiff")
target_link_libraries(bindImageIO PUBLIC bindMatrix tiff)

add_library(bindmyMath STATIC ${CMAKE_SOURCE_DIR}/src/srcMath/myMath.cpp)
target_link_libraries(bindmyMath PUBLIC bindMatrix)

add_library(bindCamera STATIC ${CMAKE_SOURCE_DIR}/src/srcMath/Camera.cpp)
target_link_libraries(bindCamera PUBLIC bindmyMath)


# Create pybind11 module
set(BINDINGS_SRC
    ${CMAKE_SOURCE_DIR}/src/bindings.cpp
    ${CMAKE_SOURCE_DIR}/src/pybind_OpenLPT/pyMatrix.cpp
    ${CMAKE_SOURCE_DIR}/src/pybind_OpenLPT/pyImageIO.cpp
    ${CMAKE_SOURCE_DIR}/src/pybind_OpenLPT/pyCamera.cpp
)
set(BINDINGS_LIB
    bindMatrix
    bindImageIO
    bindCamera
)
pybind11_add_module(pyOpenLPT ${BINDINGS_SRC})
target_link_libraries(pyOpenLPT PRIVATE ${BINDINGS_LIB})

