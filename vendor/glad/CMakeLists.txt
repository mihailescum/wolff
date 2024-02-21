# Source Files
set (SOURCES src/glad.c)
# Header Files
set (HEADER include/glad/glad.h include/KHR/khrplatform.h)

#Create Library
add_library(glad ${SOURCES} ${HEADER})

# Include Project Header Files
target_include_directories(glad PUBLIC include)