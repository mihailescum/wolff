
add_library(wolff_examples INTERFACE)

target_include_directories(wolff_examples INTERFACE include)

add_executable(simple_ising simple_ising.cpp)
target_link_libraries(simple_ising wolff)

add_subdirectory(src)

