
set (sources
     main.cpp
)
     
add_hpx_executable(hpx_bgl
SOURCES ${sources}
FOLDER
COMPONENT_DEPENDENCIES
DEPENDENCIES metis)