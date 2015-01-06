
set (sources
     headers.hpp
     main.cpp
     structures.cpp
)
     
add_executable(hpx_bgl ${sources})

target_link_libraries(hpx_bgl metis)