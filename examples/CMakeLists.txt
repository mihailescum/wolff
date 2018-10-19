
add_executable(ising ising.cpp)
add_executable(ising_animation ising_animation.cpp)
add_executable(ising_standalone ising_standalone.cpp)
add_executable(xy On.cpp)

add_compile_definitions(xy WOLFF_N=2)

target_link_libraries(ising wolff)
target_link_libraries(ising_animation wolff GL GLU glut)
target_link_libraries(ising_standalone wolff)
target_link_libraries(xy wolff)

