
add_executable(ising ising.cpp)
add_executable(ising_animation ising_animation.cpp)
add_executable(ising_standalone ising_standalone.cpp)
add_executable(xy On.cpp)
add_executable(potts_3 potts.cpp)
add_executable(clock_5 clock.cpp)
add_executable(discrete_gaussian discrete_gaussian.cpp)
add_executable(continuous_gaussian continuous_gaussian.cpp)

target_compile_definitions(xy PUBLIC WOLFF_N=2)
target_compile_definitions(potts_3 PUBLIC WOLFF_POTTSQ=3)
target_compile_definitions(clock_5 PUBLIC WOLFF_POTTSQ=5)

target_link_libraries(ising wolff)
target_link_libraries(ising_animation wolff GL GLU glut)
target_link_libraries(ising_standalone wolff)
target_link_libraries(xy wolff)
target_link_libraries(potts_3 wolff)
target_link_libraries(clock_5 wolff)
target_link_libraries(discrete_gaussian wolff)
target_link_libraries(continuous_gaussian wolff)

