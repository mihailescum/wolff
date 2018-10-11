
cmake_minimum_required(VERSION 3.9)
project(wolff DESCRIPTION "a library for preforming the wolff algorithm")

set(CMAKE_CXX_FLAGS_DEBUG "-g -Wall")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

set (CMAKE_CXX_STANDARD 17)
set (CMAKE_C_STANDARD 11)

include(GNUInstallDirs)

add_subdirectory(lib)
add_subdirectory(examples)

