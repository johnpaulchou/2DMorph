cmake_minimum_required(VERSION 3.0)
project(2DMorph)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_FIND_DEBUG_MODE 1)

# Find ROOT package
find_package(ROOT 6.30 REQUIRED COMPONENTS RooFit RooFitCore PATHS /Users/johnpaul/Work/root.6.30.06-install/cmake NO_DEFAULT_PATH)
#find_package(ROOT 6.30 REQUIRED COMPONENTS RooFit RooFitCore)
if(ROOT_FOUND)
    message(STATUS "ROOT found. Version: ${ROOT_VERSION}")
else()
    message(FATAL_ERROR "ROOT not found. Please check your ROOT installation.")
endif()


# Include ROOT headers
include_directories(${ROOT_INCLUDE_DIRS})

# Create a sources variable with a link to all cpp files to compile
set(SOURCES
    src/roofitmain.cpp
#    src/testpdfmain.cpp
    src/TH1Pdf.cc
)

# Add executable
add_executable(2DMorph ${SOURCES})

# Set the directories that should be included in the build command for this target
# when running g++ these will be included as -I/directory/path/
target_include_directories(2DMorph
    PRIVATE
        ${PROJECT_SOURCE_DIR}/include
)



# Link ROOT libraries
target_link_libraries(2DMorph PRIVATE ${ROOT_LIBRARIES})

