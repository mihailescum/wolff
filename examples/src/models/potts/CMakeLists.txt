
add_executable(wolff_3potts             wolff_potts.cpp)
add_executable(wolff_4potts             wolff_potts.cpp)
add_executable(wolff_7potts             wolff_potts.cpp)
add_executable(wolff_3clock             wolff_clock.cpp)
add_executable(wolff_5clock             wolff_clock.cpp)

set_target_properties(wolff_3potts          PROPERTIES COMPILE_FLAGS "-DPOTTSQ=3")
set_target_properties(wolff_4potts          PROPERTIES COMPILE_FLAGS "-DPOTTSQ=4")
set_target_properties(wolff_7potts          PROPERTIES COMPILE_FLAGS "-DPOTTSQ=7")
set_target_properties(wolff_3clock          PROPERTIES COMPILE_FLAGS "-DPOTTSQ=3")
set_target_properties(wolff_5clock          PROPERTIES COMPILE_FLAGS "-DPOTTSQ=5")

find_library(GL NAMES GL)
find_library(GLU NAMES GLU)
find_library(GLUT NAMES glut)

if (${GLUT} MATCHES "GLUT-NOTFOUND")
  target_link_libraries(wolff_3potts wolff wolff_examples)
  target_link_libraries(wolff_4potts wolff wolff_examples)
  target_link_libraries(wolff_7potts wolff wolff_examples)
  target_link_libraries(wolff_3clock wolff wolff_examples)
  target_link_libraries(wolff_5clock wolff wolff_examples)
else()
  target_compile_definitions(wolff_3potts PUBLIC HAVE_GLUT)
  target_compile_definitions(wolff_4potts PUBLIC HAVE_GLUT)
  target_compile_definitions(wolff_7potts PUBLIC HAVE_GLUT)
  target_compile_definitions(wolff_3clock PUBLIC HAVE_GLUT)
  target_compile_definitions(wolff_5clock PUBLIC HAVE_GLUT)

  target_link_libraries(wolff_3potts             wolff wolff_examples glut GL GLU)
  target_link_libraries(wolff_4potts             wolff wolff_examples glut GL GLU)
  target_link_libraries(wolff_7potts             wolff wolff_examples glut GL GLU)
  target_link_libraries(wolff_3clock             wolff wolff_examples glut GL GLU)
  target_link_libraries(wolff_5clock             wolff wolff_examples glut GL GLU)
endif()

install(TARGETS wolff_3potts wolff_4potts wolff_7potts wolff_3clock wolff_5clock DESTINATION ${CMAKE_INSTALL_BINDIR})

