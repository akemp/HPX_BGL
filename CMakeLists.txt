# Copyright (c) 2013 Andrew Kemp
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
# Require a recent version of cmake
cmake_minimum_required(VERSION 2.8.4 FATAL_ERROR)

# This project is C++ based.
project(hpx_BGL)

include_directories($ENV{BOOST_ROOT})
link_directories($ENV{BOOST_ROOT}/lib)

add_subdirectory(metis)
add_subdirectory(hpx_bgl)
