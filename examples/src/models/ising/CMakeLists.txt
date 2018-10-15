
add_executable(wolff_ising             wolff_ising.cpp)
add_executable(wolff_ising_2d          wolff_ising.cpp)
add_executable(wolff_ising_2d_no-field wolff_ising.cpp)
add_executable(wolff_random-field_ising wolff_random-field_ising.cpp)

set_target_properties(wolff_ising_2d          PROPERTIES COMPILE_FLAGS "-DDIMENSION=2")
set_target_properties(wolff_ising_2d_no-field PROPERTIES COMPILE_FLAGS "-DDIMENSION=2 -DNOFIELD")

find_library(GL NAMES GL)
find_library(GLU NAMES GLU)
find_library(GLUT NAMES glut)

if (${GLUT} MATCHES "GLUT-NOTFOUND")
  target_link_libraries(wolff_ising             wolff wolff_examples)
  target_link_libraries(wolff_ising_2d          wolff wolff_examples)
  target_link_libraries(wolff_ising_2d_no-field wolff wolff_examples)
  target_link_libraries(wolff_random-field_ising wolff wolff_examples)
else()
  target_compile_definitions(wolff_ising             PUBLIC HAVE_GLUT)
  target_compile_definitions(wolff_ising_2d          PUBLIC HAVE_GLUT)
  target_compile_definitions(wolff_ising_2d_no-field PUBLIC HAVE_GLUT)
  target_compile_definitions(wolff_random-field_ising PUBLIC HAVE_GLUT)

  target_link_libraries(wolff_ising             wolff wolff_examples glut GL GLU)
  target_link_libraries(wolff_ising_2d          wolff wolff_examples glut GL GLU)
  target_link_libraries(wolff_ising_2d_no-field wolff wolff_examples glut GL GLU)
  target_link_libraries(wolff_random-field_ising wolff wolff_examples glut GL GLU)
endif()

install(TARGETS wolff_ising wolff_ising_2d wolff_ising_2d_no-field wolff_random-field_ising DESTINATION ${CMAKE_INSTALL_BINDIR} OPTIONAL)

