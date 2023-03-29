cmake_minimum_required(VERSION 3.13)  # CMake version check

project(OpenLPT)               # Create project "simple_example"
set(CMAKE_CXX_STANDARD 11)         # Enable c++11 standard

# set(CMAKE_C_COMPILER "D:\\MinGW\\mingw64\\bin\\gcc.exe")
# set(CMAKE_CXX_COMPILER "D:\\MinGW\\mingw64\\bin\\g++.exe")

# # Set the output folder where your program will be created
# set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/build)
# set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/build)
# set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# # The following folder will be included
# include_directories("${PROJECT_SOURCE_DIR}")

# set(CMAKE_BUILD_TYPE "Debug")
# set(CMAKE_BUILD_TYPE "Release")

set(CMAKE_C_FLAGS_DEBUG     "${CMAKE_C_FLAGS_DEBUG}     -o0 -Wall -g -ggdb")
set(CMAKE_CXX_FLAGS_DEBUG   "${CMAKE_CXX_FLAGS_DEBUG}   -o0 -Wall -g -ggdb")
set(CMAKE_C_FLAGS_RELEASE   "${CMAKE_C_FLAGS_RELEASE}   -o3 -Wall         ")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -o3 -Wall         ")

set(HEADER_FILES
    ${CMAKE_HOME_DIRECTORY}/inc/libIO/myIO.h
    ${CMAKE_HOME_DIRECTORY}/inc/libMath/myMATH.h
    ${CMAKE_HOME_DIRECTORY}/inc/libMath/STBCommons.h
    ${CMAKE_HOME_DIRECTORY}/inc/libMath/Matrix.h
    ${CMAKE_HOME_DIRECTORY}/inc/libMath/Matrix_bool.h
    ${CMAKE_HOME_DIRECTORY}/inc/libSTB/Camera.h
    ${CMAKE_HOME_DIRECTORY}/inc/libSTB/ObjectInfo.h
    ${CMAKE_HOME_DIRECTORY}/inc/libSTB/ImageIO.h
    # ${CMAKE_HOME_DIRECTORY}/inc/libSTB/ObjectFinder.h
    # ${CMAKE_HOME_DIRECTORY}/inc/libSTB/StereoMatch.h
    # ${CMAKE_HOME_DIRECTORY}/inc/libSTB/OTF.h
    # ${CMAKE_HOME_DIRECTORY}/inc/libSTB/Shake.h
    # ${CMAKE_HOME_DIRECTORY}/inc/libSTB/IPR.h
)

set(SOURCE_FILES main.cpp
    ${CMAKE_HOME_DIRECTORY}/src/srcIO/myIO.cpp
    ${CMAKE_HOME_DIRECTORY}/src/srcMath/myMATH.cpp
    ${CMAKE_HOME_DIRECTORY}/src/srcMath/Matrix.hpp
    ${CMAKE_HOME_DIRECTORY}/src/srcMath/Matrix_bool.hpp
    ${CMAKE_HOME_DIRECTORY}/src/srcSTB/Camera.cpp
    ${CMAKE_HOME_DIRECTORY}/src/srcSTB/ObjectInfo.cpp
    ${CMAKE_HOME_DIRECTORY}/src/srcSTB/ImageIO.cpp
#     ${CMAKE_HOME_DIRECTORY}/src/srcSTB/ObjectFinder.hpp
#     ${CMAKE_HOME_DIRECTORY}/src/srcSTB/StereoMatch.hpp
#     ${CMAKE_HOME_DIRECTORY}/src/srcSTB/OTF.cpp
#     ${CMAKE_HOME_DIRECTORY}/src/srcSTB/Shake.hpp
#     ${CMAKE_HOME_DIRECTORY}/src/srcSTB/IPR.hpp
)

# Add executable target with source files listed in SOURCE_FILES variable
include_directories("${CMAKE_HOME_DIRECTORY}/inc/libIO")
include_directories("${CMAKE_HOME_DIRECTORY}/src/srcIO")
include_directories("${CMAKE_HOME_DIRECTORY}/inc/libtiff")
include_directories("${CMAKE_HOME_DIRECTORY}/inc/libMath")
include_directories("${CMAKE_HOME_DIRECTORY}/src/srcMath")
include_directories("${CMAKE_HOME_DIRECTORY}/inc/libSTB")
include_directories("${CMAKE_HOME_DIRECTORY}/src/srcSTB")

add_executable(OpenLPT ${SOURCE_FILES} ${HEADER_FILES})

add_subdirectory("${CMAKE_HOME_DIRECTORY}/inc/libtiff")
target_link_libraries(OpenLPT PUBLIC tiff)

find_package(OpenMP REQUIRED)
if (OPENMP_FOUND)
    message("OpenMP FOUND")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")

    # set_property(TARGET OpenLPT PROPERTY
    #     MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
    set_property(TARGET OpenLPT PROPERTY
        MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>DLL")
endif()
