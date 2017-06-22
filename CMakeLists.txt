
cmake_minimum_required(VERSION 3.0)
project(wolff)

include_directories(lib ~/.local/include)
link_directories(~/.local/lib)

file(GLOB SOURCES lib/*.c)
add_executable(wolff src/wolff.c ${SOURCES})

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

target_link_libraries(wolff gsl m jst cblas fftw3 tcmalloc profiler)
#target_link_libraries(pt_wolff gsl m jst cblas tcmalloc profiler)

install(TARGETS wolff DESTINATION bin)

