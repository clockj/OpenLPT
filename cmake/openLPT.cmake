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

add_library(KalmanFilter SHARED ${CMAKE_HOME_DIRECTORY}/src/srcMath/KalmanFilter.cpp)
target_link_libraries(KalmanFilter PUBLIC Matrix myMath)

add_library(ObjectInfo SHARED ${CMAKE_HOME_DIRECTORY}/src/srcObject/ObjectInfo.cpp)
target_link_libraries(ObjectInfo PUBLIC Matrix Camera)

# Sphere object 
add_library(SphereInfo SHARED ${CMAKE_HOME_DIRECTORY}/src/srcObject/Sphere/SphereInfo.cpp)
target_link_libraries(SphereInfo PUBLIC Matrix Camera)

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
target_link_libraries(Track PUBLIC Matrix myMath ObjectInfo KalmanFilter)

add_library(STB SHARED ${CMAKE_HOME_DIRECTORY}/src/srcSTB/STB.hpp)
set_target_properties(STB PROPERTIES LINKER_LANGUAGE CXX)
target_link_libraries(STB PUBLIC Matrix myMath ObjectInfo ObjectFinder StereoMatch OTF Shake IPR PredField Track KalmanFilter)

# exe
add_executable(OpenLPT src/main.cpp)
target_link_libraries(OpenLPT PRIVATE ImageIO STB Matrix KalmanFilter ObjectInfo ObjectFinder StereoMatch OTF Shake IPR PredField Track)

