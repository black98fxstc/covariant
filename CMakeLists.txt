cmake_minimum_required(VERSION 3.15)

# Ensure vcpkg is found using an absolute path relative to the project root
if(NOT DEFINED CMAKE_TOOLCHAIN_FILE)
    set(CMAKE_TOOLCHAIN_FILE "${CMAKE_CURRENT_SOURCE_DIR}/vcpkg/scripts/buildsystems/vcpkg.cmake" CACHE PATH "toolchain file")
endif()

project(CovariantProject LANGUAGES CXX)

# Set C++ standard to C++20 as required by the source code
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find FFTW3 using vcpkg configuration files. 
# Your code uses fftwf_ functions, so we specifically need the float variant.
find_package(FFTW3f CONFIG REQUIRED)

# FFTW3::fftw3f is the target provided by the vcpkg config
set(FFTW_LIBS FFTW3::fftw3f)

# The math library is needed on Unix-like systems
if(UNIX)
    find_library(MATH_LIB NAMES m REQUIRED)
endif()

# Include the current directory for headers like Weighty.hpp and cxxopts.hpp
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

# Define targets for the various tools and main programs
add_executable(testdata testdata.cpp)
target_link_libraries(testdata PRIVATE ${FFTW_LIBS} ${MATH_LIB})

add_executable(laplacian laplacian.cpp)
target_link_libraries(laplacian PRIVATE ${FFTW_LIBS} ${MATH_LIB})

add_executable(covariant_tool covariant.cpp)
target_link_libraries(covariant_tool PRIVATE ${FFTW_LIBS} ${MATH_LIB})

# Optional: Set an output name for covariant_tool to match your Makefile target
set_target_properties(covariant_tool PROPERTIES OUTPUT_NAME "covariant")