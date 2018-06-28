
cmake_minimum_required(VERSION 3.0)
project(wolff)
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

set (CMAKE_CXX_STANDARD 11)

include_directories(lib ~/.local/include)
link_directories(~/.local/lib)

file(GLOB SOURCES lib/*.c)
add_executable(wolff_finite src/wolff_finite.c ${SOURCES})
add_executable(wolff_vector src/wolff_vector.c ${SOURCES})
add_executable(wolff_dgm src/wolff_dgm.c ${SOURCES})

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

target_link_libraries(wolff_finite gsl m cblas fftw3)
target_link_libraries(wolff_vector gsl m cblas fftw3)
target_link_libraries(wolff_dgm gsl m cblas fftw3)

install(TARGETS wolff_finite wolff_vector wolff_dgm DESTINATION bin)

