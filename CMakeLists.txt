cmake_minimum_required(VERSION 3.0)
project(2DMorph)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find ROOT package
find_package(ROOT REQUIRED COMPONENTS RooFit RooFitCore)

# Include ROOT headers
include_directories(${ROOT_INCLUDE_DIRS})

# Create a sources variable with a link to all cpp files to compile
set(SOURCES
    src/main.cpp
)

# Add executable
add_executable(2DMorph src/main.cpp)

# Set the directories that should be included in the build command for this target
# when running g++ these will be included as -I/directory/path/
target_include_directories(2DMorph
    PRIVATE
        ${PROJECT_SOURCE_DIR}/include
)


# Link ROOT libraries
target_link_libraries(2DMorph ${ROOT_LIBRARIES})
