set (sources
graph_generator.cpp make_graph.cpp splittable_mrg.cpp utils.cpp
graph_generator.h make_graph.h mod_arith_32bit.h mod_arith_64bit.h
mod_arith.h mod_arith_xmt.h splittable_mrg.h utils.h user_settings.h
mrg_transitions.cpp
)
 
add_library(generator ${sources}
)