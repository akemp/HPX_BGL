# Copyright (c) 2013 Andrew Kemp
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

if(NOT HPX_ROOT)
  if ($ENV{HPX_BIN_ROOT})
    set (HPX_ROOT $ENV{HPX_BIN_ROOT})
  else()
    set (HPX_ROOT $ENV{HPX_ROOT})
  endif()
endif()
# Find and include HPX settings
if(NOT HPX_ROOT)
  message(FATAL_ERROR "HPX_ROOT not set, unable to find HPX!")
endif()
string(REPLACE "\\" "/" HPX_ROOT ${HPX_ROOT})

# Require a recent version of cmake
cmake_minimum_required(VERSION 2.8.4 FATAL_ERROR)

# This project is C++ based.
project(HPX_BGL)

# This adds the HPX cmake configuration directory to the search path.
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH}
  ${HPX_ROOT}/share/cmake-${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}/Modules)

# Instruct cmake to find the HPX settings
find_package(HPX)

include_directories(${HPX_INCLUDE_DIR})
link_directories(${HPX_LIBRARY_DIR})
include_directories($ENV{BOOST_ROOT})
link_directories($ENV{BOOST_ROOT}/lib)

# Adjust auxiliary build settings for HPX to work properly
if(NOT MSVC)
  set(HPX_RPATH ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}${HPX_RPATH})
else()
  # enable solution folders for MSVC
  set_property(GLOBAL PROPERTY USE_FOLDERS ON)
  if(HPX_NO_INSTALL)
    set(CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}/$(Configuration)" CACHE PATH
      "Prefix prepended to install directories." FORCE)
  endif()
endif()
add_subdirectory(metis)
add_subdirectory(hpx_bgl)