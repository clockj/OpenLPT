cmake_minimum_required(VERSION 3.15...3.27)  # CMake version check

project(OpenLPT)  
set(CMAKE_CXX_STANDARD 17)  


# set(CMAKE_C_COMPILER   "D:/msys2/mingw64/bin/gcc.exe")
# set(CMAKE_CXX_COMPILER "D:/msys2/mingw64/bin/g++.exe")

# set(CMAKE_BUILD_TYPE "Debug")
# set(CMAKE_BUILD_TYPE "Release")

# set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/Release)
set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)
# set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_HOME_DIRECTORY}/pyOpenLPT)
# set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_HOME_DIRECTORY}/pyOpenLPT/Release)

set(CMAKE_C_FLAGS_DEBUG     "${CMAKE_C_FLAGS_DEBUG}     -O0 -Wall -g -ggdb")
set(CMAKE_CXX_FLAGS_DEBUG   "${CMAKE_CXX_FLAGS_DEBUG}   -O0 -Wall -g -ggdb")
set(CMAKE_C_FLAGS_RELEASE   "${CMAKE_C_FLAGS_RELEASE}   -O3 -Wall         ")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -Wall         ")
# set(CMAKE_C_FLAGS_RELEASE   "${CMAKE_C_FLAGS_RELEASE}   -Ofast -Wall         ")
# set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -Ofast -Wall         ")

# Add executable target with source files listed in SOURCE_FILES variable
include_directories("${CMAKE_HOME_DIRECTORY}/inc/libtiff")
include_directories("${CMAKE_HOME_DIRECTORY}/inc/libMath")
include_directories("${CMAKE_HOME_DIRECTORY}/src/srcMath")
include_directories("${CMAKE_HOME_DIRECTORY}/inc/libObject")
include_directories("${CMAKE_HOME_DIRECTORY}/inc/libObject/Sphere")
include_directories("${CMAKE_HOME_DIRECTORY}/src/srcObject")
include_directories("${CMAKE_HOME_DIRECTORY}/src/srcObject/Sphere")
include_directories("${CMAKE_HOME_DIRECTORY}/inc/libSTB")
include_directories("${CMAKE_HOME_DIRECTORY}/src/srcSTB")

# Find openmp package
find_package(OpenMP REQUIRED)
if (OPENMP_FOUND)
    message("OpenMP FOUND")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")

    # set_property(TARGET OpenLPT PROPERTY
    #     MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
    # set_property(TARGET OpenLPT PROPERTY
    #     MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>DLL")
endif()


# Compile openLPT 
include("${CMAKE_HOME_DIRECTORY}/cmake/openLPT.cmake")

# Compile ctests 
include("${CMAKE_HOME_DIRECTORY}/cmake/ctest.cmake")


# Compile python interface
find_package(pybind11 CONFIG)
if (pybind11_FOUND)
    message("pybind11 FOUND")
    find_package(Python REQUIRED COMPONENTS Interpreter Development.Module)
    include("${CMAKE_HOME_DIRECTORY}/cmake/bindings.cmake")
else()
    message("pybind11 NOT FOUND")
endif()


# install libraries and exe
if (pybind11_FOUND)
    set(CMAKE_INSTALL_TARGETS
        Matrix
        myMath
        ImageIO
        Camera
        ObjectInfo
        SphereInfo
        ObjectFinder
        StereoMatch
        OTF
        Shake
        IPR
        PredField
        Track
        STB
        OpenLPT
        pyOpenLPT
    )
else()
    set(CMAKE_INSTALL_TARGETS
        Matrix
        myMath
        ImageIO
        Camera
        ObjectInfo
        SphereInfo
        ObjectFinder
        StereoMatch
        OTF
        Shake
        IPR
        PredField
        Track
        STB
        OpenLPT
    )
endif()

install(TARGETS ${CMAKE_INSTALL_TARGETS} DESTINATION ${CMAKE_HOME_DIRECTORY}/bin/)