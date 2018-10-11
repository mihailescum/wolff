
add_executable(wolff_dgm             wolff_dgm.cpp)
add_executable(wolff_cgm             wolff_cgm.cpp)

find_library(GL NAMES GL)
find_library(GLU NAMES GLU)
find_library(GLUT NAMES glut)

if (${GLUT} MATCHES "GLUT-NOTFOUND")
  target_link_libraries(wolff_dgm             wolff wolff_examples)
  target_link_libraries(wolff_cgm             wolff wolff_examples)
else()
  target_compile_definitions(wolff_dgm             PUBLIC HAVE_GLUT)
  target_compile_definitions(wolff_cgm             PUBLIC HAVE_GLUT)

  target_link_libraries(wolff_dgm             wolff wolff_examples glut GL GLU)
  target_link_libraries(wolff_cgm             wolff wolff_examples glut GL GLU)
endif()

install(TARGETS wolff_dgm wolff_cgm DESTINATION ${CMAKE_INSTALL_BINDIR})

