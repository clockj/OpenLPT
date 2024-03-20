cmake_minimum_required(VERSION 3.13)  # CMake version check

project(OpenLPT)               # Create project "simple_example"
set(CMAKE_CXX_STANDARD 17)         # Enable c++11 standard

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
    ${CMAKE_HOME_DIRECTORY}/inc/libSTB/Shake.h
    ${CMAKE_HOME_DIRECTORY}/inc/libSTB/IPR.h
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
    ${CMAKE_HOME_DIRECTORY}/src/srcSTB/Shake.cpp
    ${CMAKE_HOME_DIRECTORY}/src/srcSTB/IPR.cpp
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



# Build libraries
add_library(Matrix SHARED ${CMAKE_HOME_DIRECTORY}/src/srcMath/Matrix.hpp)
set_target_properties(Matrix PROPERTIES LINKER_LANGUAGE CXX)

add_library(myMath SHARED ${CMAKE_HOME_DIRECTORY}/src/srcMath/myMATH.cpp)
target_link_libraries(myMath PUBLIC Matrix)

add_library(ImageIO SHARED ${CMAKE_HOME_DIRECTORY}/src/srcMath/ImageIO.cpp)
add_subdirectory("${CMAKE_HOME_DIRECTORY}/inc/libtiff")
target_link_libraries(ImageIO PUBLIC Matrix tiff)

add_library(Camera SHARED ${CMAKE_HOME_DIRECTORY}/src/srcMath/Camera.cpp)
target_link_libraries(Camera PUBLIC Matrix myMath)

add_library(ObjectInfo SHARED ${CMAKE_HOME_DIRECTORY}/src/srcObject/ObjectInfo.cpp)
target_link_libraries(ObjectInfo PUBLIC Matrix Camera)

add_library(ObjectFinder SHARED ${CMAKE_HOME_DIRECTORY}/src/srcObject/ObjectFinder.hpp)
set_target_properties(ObjectFinder PROPERTIES LINKER_LANGUAGE CXX)
target_link_libraries(ObjectFinder PUBLIC Matrix myMath ObjectInfo)

add_library(StereoMatch SHARED ${CMAKE_HOME_DIRECTORY}/src/srcSTB/StereoMatch.hpp)
set_target_properties(StereoMatch PROPERTIES LINKER_LANGUAGE CXX)
target_link_libraries(StereoMatch PUBLIC Matrix myMath ObjectInfo Camera)

add_library(OTF SHARED ${CMAKE_HOME_DIRECTORY}/src/srcSTB/OTF.cpp)
target_link_libraries(OTF PUBLIC Matrix myMath)

add_library(Shake SHARED ${CMAKE_HOME_DIRECTORY}/src/srcSTB/Shake.cpp)
target_link_libraries(Shake PUBLIC Matrix myMath ObjectInfo Camera OTF)

add_library(IPR SHARED ${CMAKE_HOME_DIRECTORY}/src/srcSTB/IPR.hpp)
set_target_properties(IPR PROPERTIES LINKER_LANGUAGE CXX)
target_link_libraries(IPR PUBLIC Matrix Camera ObjectInfo ObjectFinder StereoMatch Shake OTF)

add_library(PredField SHARED ${CMAKE_HOME_DIRECTORY}/src/srcSTB/PredField.hpp)
set_target_properties(PredField PROPERTIES LINKER_LANGUAGE CXX)
target_link_libraries(PredField PUBLIC Matrix myMath ObjectInfo)

add_library(Track SHARED ${CMAKE_HOME_DIRECTORY}/src/srcSTB/Track.hpp)
set_target_properties(Track PROPERTIES LINKER_LANGUAGE CXX)
target_link_libraries(Track PUBLIC Matrix myMath ObjectInfo)

add_library(STB SHARED ${CMAKE_HOME_DIRECTORY}/src/srcSTB/STB.hpp)
set_target_properties(STB PROPERTIES LINKER_LANGUAGE CXX)
target_link_libraries(STB PUBLIC Matrix myMath ObjectInfo ObjectFinder StereoMatch OTF Shake IPR PredField Track)


# install libraries
set(CMAKE_INSTALL_TARGETS
    Matrix
    myMath
    ImageIO
    Camera
    ObjectInfo
    ObjectFinder
    StereoMatch
    OTF
    Shake
    IPR
    PredField
    Track
    STB
)
install(TARGETS ${CMAKE_INSTALL_TARGETS} DESTINATION ${CMAKE_HOME_DIRECTORY}/bin/)

# # add_executable(OpenLPT ${SOURCEFILES} ${HEADERFILES})
# add_executable(OpenLPT main.cpp)
# target_link_libraries(OpenLPT PUBLIC libMatrix)

# add_subdirectory("${CMAKE_HOME_DIRECTORY}/inc/libtiff")
# target_link_libraries(OpenLPT PUBLIC tiff)


# Enable testing with CTest
enable_testing()

# test header
set(TESTHEADERFILES
    ${CMAKE_HOME_DIRECTORY}/test/test.h
)

# test Matrix
add_executable(test_Matrix 
    ${CMAKE_HOME_DIRECTORY}/test/test_Matrix.cpp
    ${TESTHEADERFILES}
)
target_link_libraries(test_Matrix PRIVATE Matrix)
add_test(test_Matrix test_Matrix)
set_tests_properties(test_Matrix
                     PROPERTIES FAIL_REGULAR_EXPRESSION "failed on line")


# test myMATH
add_executable(test_myMATH 
    ${CMAKE_HOME_DIRECTORY}/test/test_myMATH.cpp
    ${TESTHEADERFILES}
)
target_link_libraries(test_myMATH PRIVATE myMath Matrix)
add_test(test_myMATH test_myMATH)
set_tests_properties(test_myMATH
                     PROPERTIES FAIL_REGULAR_EXPRESSION "failed on line")


# test ImageIO
add_executable(test_ImageIO 
    ${CMAKE_HOME_DIRECTORY}/test/test_ImageIO.cpp
    ${TESTHEADERFILES}
)
target_link_libraries(test_ImageIO PRIVATE ImageIO Matrix tiff)
add_test(test_ImageIO test_ImageIO)
set_tests_properties(test_ImageIO
                     PROPERTIES FAIL_REGULAR_EXPRESSION "failed on line")


# test Camera
add_executable(test_Camera 
    ${CMAKE_HOME_DIRECTORY}/test/test_Camera.cpp
    ${TESTHEADERFILES}
)
target_link_libraries(test_Camera PRIVATE Camera Matrix myMath)
add_test(test_Camera test_Camera)
set_tests_properties(test_Camera
                     PROPERTIES FAIL_REGULAR_EXPRESSION "failed on line")


# test ObjectInfo
add_executable(test_ObjectInfo 
    ${CMAKE_HOME_DIRECTORY}/test/test_ObjectInfo.cpp
    ${TESTHEADERFILES}
)
target_link_libraries(test_ObjectInfo PRIVATE ObjectInfo Camera)
add_test(test_ObjectInfo test_ObjectInfo)
set_tests_properties(test_ObjectInfo
                     PROPERTIES FAIL_REGULAR_EXPRESSION "failed on line")


# test ObjectFinder
add_executable(test_ObjectFinder 
    ${CMAKE_HOME_DIRECTORY}/test/test_ObjectFinder.cpp
    ${TESTHEADERFILES}
)
target_link_libraries(test_ObjectFinder PRIVATE ObjectFinder ObjectInfo ImageIO)
add_test(test_ObjectFinder test_ObjectFinder)
set_tests_properties(test_ObjectFinder
                     PROPERTIES FAIL_REGULAR_EXPRESSION "failed on line")


# test StereoMatch
add_executable(test_StereoMatch 
    ${CMAKE_HOME_DIRECTORY}/test/test_StereoMatch.cpp
    ${TESTHEADERFILES}
)
target_link_libraries(test_StereoMatch PRIVATE Camera StereoMatch ObjectInfo ImageIO ObjectFinder myMath Matrix)
add_test(test_StereoMatch test_StereoMatch)
set_tests_properties(test_StereoMatch
                     PROPERTIES FAIL_REGULAR_EXPRESSION "failed on line")


# test OTF
add_executable(test_OTF 
    ${CMAKE_HOME_DIRECTORY}/test/test_OTF.cpp
    ${TESTHEADERFILES}
)
target_link_libraries(test_OTF PRIVATE OTF Matrix)
add_test(test_OTF test_OTF)
set_tests_properties(test_OTF
                     PROPERTIES FAIL_REGULAR_EXPRESSION "failed on line")


# test Shake
add_executable(test_Shake 
    ${CMAKE_HOME_DIRECTORY}/test/test_Shake.cpp
    ${TESTHEADERFILES}
)
target_link_libraries(test_Shake PRIVATE Shake Shake ObjectFinder StereoMatch ObjectInfo Camera ImageIO OTF myMath Matrix)
add_test(test_Shake test_Shake)
set_tests_properties(test_Shake
                     PROPERTIES FAIL_REGULAR_EXPRESSION "failed on line")


# test IPR
add_executable(test_IPR 
    ${CMAKE_HOME_DIRECTORY}/test/test_IPR.cpp
    ${TESTHEADERFILES}
)
target_link_libraries(test_IPR PRIVATE 
    Matrix
    myMath
    Camera
    ImageIO
    ObjectInfo
    ObjectFinder
    StereoMatch
    Shake
    OTF
    IPR
)
add_test(test_IPR test_IPR)
set_tests_properties(test_IPR
                     PROPERTIES FAIL_REGULAR_EXPRESSION "failed on line")


# test PredField
add_executable(test_PredField 
    ${CMAKE_HOME_DIRECTORY}/test/test_PredField.cpp
    ${TESTHEADERFILES}
)
target_link_libraries(test_PredField PRIVATE PredField Matrix ObjectInfo)
add_test(test_PredField test_PredField)


# test Track
add_executable(test_Track 
    ${CMAKE_HOME_DIRECTORY}/test/test_Track.cpp
    ${TESTHEADERFILES}
)
target_link_libraries(test_Track PRIVATE Track Matrix ObjectInfo)
add_test(test_Track test_Track)


# test STB
add_executable(test_STB 
    ${CMAKE_HOME_DIRECTORY}/test/test_STB.cpp
    ${TESTHEADERFILES}
)
target_link_libraries(test_STB PRIVATE ImageIO STB Matrix ObjectInfo ObjectFinder StereoMatch OTF Shake IPR PredField Track)
add_test(test_STB test_STB)