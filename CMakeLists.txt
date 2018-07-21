
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

find_library(GSL gsl)
find_library(FFTW fftw3)
find_library(M m)
FIND_LIBRARY(GL GL)
FIND_LIBRARY(GLU GLU)
FIND_LIBRARY(GLUT glut)

target_link_libraries(wolff_finite cblas gsl m)
target_link_libraries(analyze_correlations cblas gsl fftw3 m)
if (GL_FOUND AND GLU_FOUND AND GLUT_FOUND)
  target_link_libraries(wolff_ising cblas gsl m glut GL GLU)
  target_link_libraries(wolff_heisenberg cblas gsl m glut GL GLU)
  target_link_libraries(wolff_planar cblas gsl m glut GL GLU)
  target_compile_definitions(wolff_ising HAVE_GLUT)
  target_compile_definitions(wolff_planar HAVE_GLUT)
  target_compile_definitions(wolff_heisenberg HAVE_GLUT)
else ()
  target_link_libraries(wolff_ising cblas gsl m)
  target_link_libraries(wolff_heisenberg cblas gsl m)
  target_link_libraries(wolff_planar cblas gsl m)
endif()

install(TARGETS wolff_finite wolff_ising wolff_heisenberg wolff_planar analyze_correlations DESTINATION bin)

