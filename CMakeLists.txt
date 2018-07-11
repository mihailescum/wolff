
cmake_minimum_required(VERSION 3.0)
project(wolff)
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

set (CMAKE_CXX_STANDARD 17)

include_directories(lib ~/.local/include)
link_directories(~/.local/lib)

file(GLOB CSOURCES lib/*.c)
file(GLOB CPPSOURCES lib/*.cpp)

add_executable(wolff_finite src/wolff_finite.c ${CSOURCES})
add_executable(wolff_heisenberg src/wolff_heisenberg.cpp ${CPPSOURCES} ${CSOURCES})
add_executable(wolff_planar src/wolff_planar.cpp ${CPPSOURCES} ${CSOURCES})
add_executable(analyze_correlations src/analyze_correlations.cpp ${CPPSOURCES} ${CSOURCES})

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

target_link_libraries(wolff_finite gsl m cblas fftw3 tcmalloc profiler)
target_link_libraries(wolff_heisenberg gsl m cblas fftw3 tcmalloc profiler)
target_link_libraries(wolff_planar gsl m cblas fftw3 tcmalloc profiler)
target_link_libraries(analyze_correlations gsl m cblas fftw3 profiler)

install(TARGETS wolff_finite wolff_heisenberg wolff_planar analyze_correlations DESTINATION bin)

