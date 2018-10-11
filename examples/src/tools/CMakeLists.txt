
find_library(fftw REQUIRED NAMES fftw3)

add_executable(analyze_correlations analyze_correlations.cpp)

target_link_libraries(analyze_correlations fftw3 wolff)

install(TARGETS analyze_correlations DESTINATION ${CMAKE_INSTALL_BINDIR} OPTIONAL)

