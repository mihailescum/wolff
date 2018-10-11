
add_executable(wolff_planar               wolff_On.cpp)
add_executable(wolff_planar_2d_no-field   wolff_On.cpp)
add_executable(wolff_heisenberg           wolff_On.cpp)

set_target_properties(wolff_planar             PROPERTIES COMPILE_FLAGS "-DN_COMP=2")
set_target_properties(wolff_planar_2d_no-field PROPERTIES COMPILE_FLAGS "-DN_COMP=2 -DDIMENSION=2 -DNOFIELD")
set_target_properties(wolff_heisenberg         PROPERTIES COMPILE_FLAGS "-DN_COMP=3")

find_library(GL NAMES GL)
find_library(GLU NAMES GLU)
find_library(GLUT NAMES glut)

if (${GLUT} MATCHES "GLUT-NOTFOUND")
  target_link_libraries(wolff_planar             wolff wolff_examples)
  target_link_libraries(wolff_planar_2d_no-field wolff wolff_examples)
  target_link_libraries(wolff_heisenberg         wolff wolff_examples)
else()
  target_compile_definitions(wolff_planar             PUBLIC HAVE_GLUT)
  target_compile_definitions(wolff_planar_2d_no-field             PUBLIC HAVE_GLUT)
  target_compile_definitions(wolff_heisenberg             PUBLIC HAVE_GLUT)

  target_link_libraries(wolff_planar             wolff wolff_examples glut GL GLU)
  target_link_libraries(wolff_planar_2d_no-field wolff wolff_examples glut GL GLU)
  target_link_libraries(wolff_heisenberg             wolff wolff_examples glut GL GLU)
endif()

install(TARGETS wolff_planar wolff_planar_2d_no-field wolff_heisenberg DESTINATION ${CMAKE_INSTALL_BINDIR} OPTIONAL)

