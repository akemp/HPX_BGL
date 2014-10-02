
set (sources
     main.cpp
     validate.cpp
)
     
add_hpx_executable(hpx_bgl
SOURCES ${sources}
FOLDER
COMPONENT_DEPENDENCIES
DEPENDENCIES graph500)