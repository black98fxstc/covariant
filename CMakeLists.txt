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
# vcpkg's fftw3 port provides separate config files for float (FFTW3f) and double (FFTW3) precisions.
# The code primarily uses fftwf_ functions, but loading both the base (double) and float 
# configurations ensures that all threaded targets in the FFTW3:: namespace are found.
find_package(FFTW3 CONFIG REQUIRED)
find_package(FFTW3f CONFIG REQUIRED)

find_package(Matplot++ CONFIG REQUIRED)
find_package(cxxopts CONFIG REQUIRED)
find_package(nlohmann_json CONFIG REQUIRED)
# Find the packages provided by vcpkg
find_package(LibXml2 REQUIRED)
find_package(LibXslt REQUIRED)

# Note: vcpkg provides FFTW3::fftw3f and FFTW3::fftw3f_threads when the threads feature is enabled.
# We will use these directly in target_link_libraries
 
# The math library is needed on Unix-like systems
if(UNIX)
    find_library(MATH_LIB NAMES m REQUIRED)
endif()

# Include the current directory for headers like Weighty.hpp and cxxopts.hpp
include_directories(${CMAKE_CURRENT_SOURCE_DIR} /usr/include/libxml2/ )

# Define targets for the various tools and main programs
add_executable(testdata testdata.cpp)
# testdata is a generator and does not use the Weighty class's FFTW-based filtering.
target_link_libraries(testdata PRIVATE cxxopts::cxxopts ${MATH_LIB})

add_executable(laplace laplace.cpp)
target_link_libraries(laplace PRIVATE FFTW3::fftw3f cxxopts::cxxopts ${MATH_LIB})

add_executable(riemann riemann.cpp)
target_link_libraries(riemann PRIVATE FFTW3::fftw3f cxxopts::cxxopts ${MATH_LIB})

add_executable(darwin darwin.cpp)
target_link_libraries(darwin PRIVATE FFTW3::fftw3f Matplot++::matplot nlohmann_json::nlohmann_json cxxopts::cxxopts ${MATH_LIB} LibXml2::LibXml2 LibXslt::LibXslt )
# Add LibXml2::LibXml2 and LibXslt::LibXslt to your existing link libraries

# Set output names to match your Makefile targets
set_target_properties(riemann PROPERTIES OUTPUT_NAME "riemann")
set_target_properties(darwin PROPERTIES OUTPUT_NAME "darwin")