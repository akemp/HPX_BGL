# Copyright (c) 2011 Matt Anderson
# Copyright (c) 2007-2012 Hartmut Kaiser
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

set (sources
 graph500.c   kronecker.h        options.h  rmat.cpp  xalloc.c
 compat.h        graph500.h make-edgelist.h  make-edgelist.cpp  prng.cpp   rmat.h
 kronecker.c          prng.h     verify.c verify.h
 generator/utils.h generator/utils.c generator/splittable_mrg.h
 generator/splittable_mrg.c generator/make_graph.h generator/make_graph.c
 generator/graph_generator.h generator/graph_generator.c
)
     
add_library(graph500 ${sources})
include_directories (${CMAKE_CURRENT_SOURCE_DIR})