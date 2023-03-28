# Install script for directory: D:/JHU Codes/OpenLPT/OpenLPT/inc/libtiff

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "C:/Program Files (x86)/OpenLPT")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Debug")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

# Set default install directory permissions.
if(NOT DEFINED CMAKE_OBJDUMP)
  set(CMAKE_OBJDUMP "D:/MinGW/mingw64/bin/objdump.exe")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "D:/JHU Codes/OpenLPT/OpenLPT/inc/libtiff/install/lib/libtiff.a")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "D:/JHU Codes/OpenLPT/OpenLPT/inc/libtiff/install/lib" TYPE STATIC_LIBRARY FILES "D:/JHU Codes/OpenLPT/OpenLPT/build/inc/libtiff/libtiff.a")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "D:/JHU Codes/OpenLPT/OpenLPT/inc/libtiff/install/include/tiff.h;D:/JHU Codes/OpenLPT/OpenLPT/inc/libtiff/install/include/tiffio.h;D:/JHU Codes/OpenLPT/OpenLPT/inc/libtiff/install/include/tiffvers.h;D:/JHU Codes/OpenLPT/OpenLPT/inc/libtiff/install/include/tiffconf.h")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "D:/JHU Codes/OpenLPT/OpenLPT/inc/libtiff/install/include" TYPE FILE FILES
    "D:/JHU Codes/OpenLPT/OpenLPT/inc/libtiff/tiff.h"
    "D:/JHU Codes/OpenLPT/OpenLPT/inc/libtiff/tiffio.h"
    "D:/JHU Codes/OpenLPT/OpenLPT/inc/libtiff/tiffvers.h"
    "D:/JHU Codes/OpenLPT/OpenLPT/build/inc/libtiff/tiffconf.h"
    )
endif()

