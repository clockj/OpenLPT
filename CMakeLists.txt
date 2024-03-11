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

set(CMAKE_C_FLAGS_DEBUG     "${CMAKE_C_FLAGS_DEBUG}     -O0 -Wall -g -ggdb")
set(CMAKE_CXX_FLAGS_DEBUG   "${CMAKE_CXX_FLAGS_DEBUG}   -O0 -Wall -g -ggdb")
set(CMAKE_C_FLAGS_RELEASE   "${CMAKE_C_FLAGS_RELEASE}   -O3 -Wall         ")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -Wall         ")
# set(CMAKE_C_FLAGS_RELEASE   "${CMAKE_C_FLAGS_RELEASE}   -Ofast -Wall         ")
# set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -Ofast -Wall         ")

set(HEADERFILES
    ${CMAKE_HOME_DIRECTORY}/inc/libMath/myMATH.h
    ${CMAKE_HOME_DIRECTORY}/inc/libMath/STBCommons.h
    ${CMAKE_HOME_DIRECTORY}/inc/libMath/Matrix.h
    ${CMAKE_HOME_DIRECTORY}/inc/libMath/ImageIO.h
    ${CMAKE_HOME_DIRECTORY}/inc/libMath/Camera.h
    ${CMAKE_HOME_DIRECTORY}/inc/libObject/ObjectInfo.h
    ${CMAKE_HOME_DIRECTORY}/inc/libObject/ObjectFinder.h
    ${CMAKE_HOME_DIRECTORY}/inc/libSTB/StereoMatch.h
    ${CMAKE_HOME_DIRECTORY}/inc/libSTB/OTF.h
    # ${CMAKE_HOME_DIRECTORY}/inc/libSTB/Shake.h
    # ${CMAKE_HOME_DIRECTORY}/inc/libSTB/IPR.h
    # ${CMAKE_HOME_DIRECTORY}/inc/libSTB/PredField.h
    # ${CMAKE_HOME_DIRECTORY}/inc/libSTB/Track.h
    # ${CMAKE_HOME_DIRECTORY}/inc/libSTB/STB.h
)

set(SOURCEFILES main.cpp
    ${CMAKE_HOME_DIRECTORY}/src/srcMath/myMATH.cpp
    ${CMAKE_HOME_DIRECTORY}/src/srcMath/Matrix.hpp
    ${CMAKE_HOME_DIRECTORY}/src/srcMath/ImageIO.cpp
    ${CMAKE_HOME_DIRECTORY}/src/srcMath/Camera.cpp
    ${CMAKE_HOME_DIRECTORY}/src/srcObject/ObjectInfo.cpp
    ${CMAKE_HOME_DIRECTORY}/src/srcObject/ObjectFinder.hpp
    ${CMAKE_HOME_DIRECTORY}/src/srcSTB/StereoMatch.hpp
    ${CMAKE_HOME_DIRECTORY}/src/srcSTB/OTF.cpp
    # ${CMAKE_HOME_DIRECTORY}/src/srcSTB/Shake.hpp
    # ${CMAKE_HOME_DIRECTORY}/src/srcSTB/IPR.hpp
    # ${CMAKE_HOME_DIRECTORY}/src/srcSTB/PredField.hpp
    # ${CMAKE_HOME_DIRECTORY}/src/srcSTB/Track.hpp
    # ${CMAKE_HOME_DIRECTORY}/src/srcSTB/STB.hpp
)

# Add executable target with source files listed in SOURCE_FILES variable
include_directories("${CMAKE_HOME_DIRECTORY}/inc/libtiff")
include_directories("${CMAKE_HOME_DIRECTORY}/inc/libMath")
include_directories("${CMAKE_HOME_DIRECTORY}/src/srcMath")
include_directories("${CMAKE_HOME_DIRECTORY}/inc/libObject")
include_directories("${CMAKE_HOME_DIRECTORY}/src/srcObject")
include_directories("${CMAKE_HOME_DIRECTORY}/inc/libSTB")
include_directories("${CMAKE_HOME_DIRECTORY}/src/srcSTB")

add_executable(OpenLPT ${SOURCEFILES} ${HEADERFILES})

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



# Enable testing with CTest
enable_testing()

# test
set(TESTHEADERFILES
    ${CMAKE_HOME_DIRECTORY}/test/test.h
)


# test Matrix
add_executable(test_Matrix ${CMAKE_HOME_DIRECTORY}/test/test_Matrix.cpp
    ${CMAKE_HOME_DIRECTORY}/src/srcMath/Matrix.hpp
    ${TESTHEADERFILES}
    ${CMAKE_HOME_DIRECTORY}/inc/libMath/Matrix.h
    ${CMAKE_HOME_DIRECTORY}/inc/libMath/STBCommons.h
)
add_test(test_Matrix test_Matrix)
set_tests_properties(test_Matrix
                     PROPERTIES FAIL_REGULAR_EXPRESSION "failed on line")
# add_test(test_Matrix_compare ${CMAKE_COMMAND} -E compare_files ${CMAKE_HOME_DIRECTORY}/test/results/test_Matrix/test_function_5.csv ${CMAKE_HOME_DIRECTORY}/test/solutions/test_Matrix/test_function_5.csv)


# test myMATH
add_executable(test_myMATH ${CMAKE_HOME_DIRECTORY}/test/test_myMATH.cpp
    ${CMAKE_HOME_DIRECTORY}/src/srcMath/myMATH.cpp
    ${TESTHEADERFILES}
    ${CMAKE_HOME_DIRECTORY}/inc/libMath/myMATH.h
)
add_test(test_myMATH test_myMATH)
set_tests_properties(test_myMATH
                     PROPERTIES FAIL_REGULAR_EXPRESSION "failed on line")


# test ImageIO
add_executable(test_ImageIO ${CMAKE_HOME_DIRECTORY}/test/test_ImageIO.cpp
    ${CMAKE_HOME_DIRECTORY}/src/srcMath/ImageIO.cpp
    ${TESTHEADERFILES}
    ${CMAKE_HOME_DIRECTORY}/inc/libMath/ImageIO.h
)
target_link_libraries(test_ImageIO PUBLIC tiff)
add_test(test_ImageIO test_ImageIO)
set_tests_properties(test_ImageIO
                     PROPERTIES FAIL_REGULAR_EXPRESSION "failed on line")


# test Camera
add_executable(test_Camera ${CMAKE_HOME_DIRECTORY}/test/test_Camera.cpp
    ${CMAKE_HOME_DIRECTORY}/src/srcMath/Camera.cpp
    ${CMAKE_HOME_DIRECTORY}/src/srcMath/myMATH.cpp
    ${TESTHEADERFILES}
    ${CMAKE_HOME_DIRECTORY}/inc/libMath/Camera.h
    ${CMAKE_HOME_DIRECTORY}/inc/libMath/myMATH.h
)
add_test(test_Camera test_Camera)
set_tests_properties(test_Camera
                     PROPERTIES FAIL_REGULAR_EXPRESSION "failed on line")


# test ObjectInfo
add_executable(test_ObjectInfo ${CMAKE_HOME_DIRECTORY}/test/test_ObjectInfo.cpp
    ${CMAKE_HOME_DIRECTORY}/src/srcObject/ObjectInfo.cpp
    ${CMAKE_HOME_DIRECTORY}/src/srcMath/Camera.cpp
    ${CMAKE_HOME_DIRECTORY}/src/srcMath/myMATH.cpp
    ${TESTHEADERFILES}
    ${CMAKE_HOME_DIRECTORY}/inc/libObject/ObjectInfo.h
    ${CMAKE_HOME_DIRECTORY}/inc/libMath/Camera.h
    ${CMAKE_HOME_DIRECTORY}/inc/libMath/myMATH.h
)
add_test(test_ObjectInfo test_ObjectInfo)
set_tests_properties(test_ObjectInfo
                     PROPERTIES FAIL_REGULAR_EXPRESSION "failed on line")


# test ObjectFinder
add_executable(test_ObjectFinder ${CMAKE_HOME_DIRECTORY}/test/test_ObjectFinder.cpp
    ${CMAKE_HOME_DIRECTORY}/src/srcObject/ObjectFinder.hpp
    ${CMAKE_HOME_DIRECTORY}/src/srcObject/ObjectInfo.cpp
    ${CMAKE_HOME_DIRECTORY}/src/srcMath/Camera.cpp
    ${CMAKE_HOME_DIRECTORY}/src/srcMath/myMATH.cpp
    ${CMAKE_HOME_DIRECTORY}/src/srcMath/ImageIO.cpp
    ${TESTHEADERFILES}
    ${CMAKE_HOME_DIRECTORY}/inc/libObject/ObjectFinder.h
    ${CMAKE_HOME_DIRECTORY}/inc/libObject/ObjectInfo.h
    ${CMAKE_HOME_DIRECTORY}/inc/libMath/Camera.h
    ${CMAKE_HOME_DIRECTORY}/inc/libMath/myMATH.h
    ${CMAKE_HOME_DIRECTORY}/inc/libMath/ImageIO.h
)
target_link_libraries(test_ObjectFinder PUBLIC tiff)
add_test(test_ObjectFinder test_ObjectFinder)
set_tests_properties(test_ObjectFinder
                     PROPERTIES FAIL_REGULAR_EXPRESSION "failed on line")


# test StereoMatch
add_executable(test_StereoMatch ${CMAKE_HOME_DIRECTORY}/test/test_StereoMatch.cpp
    ${CMAKE_HOME_DIRECTORY}/src/srcSTB/StereoMatch.hpp
    ${CMAKE_HOME_DIRECTORY}/src/srcMath/Matrix.hpp
    ${CMAKE_HOME_DIRECTORY}/src/srcMath/ImageIO.cpp
    ${CMAKE_HOME_DIRECTORY}/src/srcMath/Camera.cpp
    ${CMAKE_HOME_DIRECTORY}/src/srcMath/myMATH.cpp
    ${CMAKE_HOME_DIRECTORY}/src/srcObject/ObjectInfo.cpp
    ${CMAKE_HOME_DIRECTORY}/src/srcObject/ObjectFinder.hpp
    ${TESTHEADERFILES}
    ${CMAKE_HOME_DIRECTORY}/inc/libSTB/StereoMatch.h
    ${CMAKE_HOME_DIRECTORY}/inc/libMath/Matrix.h
    ${CMAKE_HOME_DIRECTORY}/inc/libMath/ImageIO.h
    ${CMAKE_HOME_DIRECTORY}/inc/libMath/Camera.h
    ${CMAKE_HOME_DIRECTORY}/inc/libMath/myMATH.h
    ${CMAKE_HOME_DIRECTORY}/inc/libObject/ObjectInfo.h
    ${CMAKE_HOME_DIRECTORY}/inc/libObject/ObjectFinder.h
)
target_link_libraries(test_StereoMatch PUBLIC tiff)
add_test(test_StereoMatch test_StereoMatch)
set_tests_properties(test_StereoMatch
                     PROPERTIES FAIL_REGULAR_EXPRESSION "failed on line")


# test OTF
add_executable(test_OTF ${CMAKE_HOME_DIRECTORY}/test/test_OTF.cpp
    ${CMAKE_HOME_DIRECTORY}/src/srcSTB/OTF.cpp
    ${CMAKE_HOME_DIRECTORY}/src/srcMath/Matrix.hpp
    ${CMAKE_HOME_DIRECTORY}/src/srcMath/myMATH.cpp
    ${TESTHEADERFILES}
    ${CMAKE_HOME_DIRECTORY}/inc/libSTB/OTF.h
    ${CMAKE_HOME_DIRECTORY}/inc/libMath/Matrix.h
    ${CMAKE_HOME_DIRECTORY}/inc/libMath/myMATH.h
)
add_test(test_OTF test_OTF)
set_tests_properties(test_OTF
                     PROPERTIES FAIL_REGULAR_EXPRESSION "failed on line")