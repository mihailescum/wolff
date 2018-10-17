
add_executable(ising ising.cpp)
add_executable(animate_ising animate_ising.cpp)
add_executable(xy On.cpp)

add_compile_definitions(xy WOLFF_N=2)

target_link_libraries(ising wolff)
target_link_libraries(animate_ising wolff GL GLU glut)
target_link_libraries(xy wolff)

