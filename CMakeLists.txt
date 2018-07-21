
cmake_minimum_required(VERSION 3.0)
project(wolff)
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

set (CMAKE_CXX_STANDARD 17)
set (CMAKE_C_STANDARD 11)

include_directories(lib ~/.local/include)
link_directories(~/.local/lib)

file(GLOB CSOURCES lib/*.c)
file(GLOB CPPSOURCES lib/*.cpp)

add_executable(wolff_finite src/wolff_finite.c ${CSOURCES})
add_executable(wolff_ising src/wolff_ising.cpp ${CPPSOURCES} ${CSOURCES})
add_executable(wolff_planar src/wolff_On.cpp ${CPPSOURCES} ${CSOURCES})
add_executable(wolff_heisenberg src/wolff_On.cpp ${CPPSOURCES} ${CSOURCES})
add_executable(analyze_correlations src/analyze_correlations.cpp ${CPPSOURCES} ${CSOURCES})

SET_TARGET_PROPERTIES(wolff_planar PROPERTIES COMPILE_FLAGS "-DN_COMP=2")
SET_TARGET_PROPERTIES(wolff_heisenberg PROPERTIES COMPILE_FLAGS "-DN_COMP=3")

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

target_link_libraries(wolff_finite gsl cblas fftw3 m)
target_link_libraries(wolff_ising gsl cblas fftw3 m glut GL GLU)
target_link_libraries(wolff_heisenberg gsl cblas fftw3 m)
target_link_libraries(wolff_planar gsl cblas fftw3 m)
target_link_libraries(analyze_correlations gsl cblas fftw3 m)

install(TARGETS wolff_finite wolff_ising wolff_heisenberg wolff_planar analyze_correlations DESTINATION bin)

