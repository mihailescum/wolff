
project(libwolff LANGUAGES CXX)

add_library(wolff INTERFACE)
target_include_directories(wolff INTERFACE .)
target_compile_features(wolff INTERFACE cxx_std_17)

install(FILES "wolff.hpp" DESTINATION include)
install(DIRECTORY "wolff_models" DESTINATION include)

