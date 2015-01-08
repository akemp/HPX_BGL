
set (sources
     headers.hpp
     main.hpp
     main.cpp
     structures.cpp
)
     
add_hpx_executable(hpx_bgl
SOURCES ${sources}
FOLDER
COMPONENT_DEPENDENCIES iostreams
DEPENDENCIES metis)