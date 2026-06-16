cmake_minimum_required(VERSION 3.15)

# Ensure vcpkg is found using an absolute path relative to the project root
if(NOT DEFINED CMAKE_TOOLCHAIN_FILE)
    set(CMAKE_TOOLCHAIN_FILE "${CMAKE_CURRENT_SOURCE_DIR}/vcpkg/scripts/buildsystems/vcpkg.cmake" CACHE PATH "toolchain file")
endif()

project(CovariantProject LANGUAGES CXX C)

# Set C++ standard to C++20 as required by the source code
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(MSVC)
    # MSVC needs /bigobj for heavily templated code like Eigen
    add_compile_options(/bigobj)
endif()

find_package(FFTW3 CONFIG REQUIRED)
find_package(FFTW3f CONFIG REQUIRED)
find_package(Matplot++ CONFIG REQUIRED)
find_package(cxxopts CONFIG REQUIRED)
find_package(nlohmann_json CONFIG REQUIRED)
find_package(Eigen3 CONFIG REQUIRED)
find_package(LibXml2 REQUIRED)
find_package(LibXslt REQUIRED)
find_package(ftxui CONFIG REQUIRED)

# Threads are required for Workers.cpp (std::thread, std::mutex, std::condition_variable)
set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)

# Note: vcpkg provides FFTW3::fftw3f and FFTW3::fftw3f_threads when the threads feature is enabled.
# We will use these directly in target_link_libraries
 
# The math library is needed on Unix-like systems
if(UNIX)
    find_library(MATH_LIB NAMES m REQUIRED)
endif()

# Include the current directory for headers like Weighty.hpp and cxxopts.hpp
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

# Create a static library for our separated compilation units so we don't compile them repeatedly
add_library(covariant_core STATIC Samples.cpp Samples.hpp Workers.cpp Workers.hpp Gating.cpp Gating.hpp Transforms.cpp Transforms.hpp Weighty.hpp)
target_link_libraries(covariant_core PUBLIC LibXml2::LibXml2 LibXslt::LibXslt Threads::Threads)

# Define targets for the various tools and main programs
add_executable(testdata testdata.cpp)
# testdata is a generator and does not use the Weighty class's FFTW-based filtering.
target_link_libraries(testdata PRIVATE covariant_core cxxopts::cxxopts ${MATH_LIB})

# add_executable(laplace laplace.cpp)
# target_link_libraries(laplace PRIVATE covariant_core FFTW3::fftw3f cxxopts::cxxopts ${MATH_LIB})

# add_executable(riemann riemann.cpp)
# target_link_libraries(riemann PRIVATE covariant_core FFTW3::fftw3f cxxopts::cxxopts ${MATH_LIB})

# add_executable(darwin darwin.cpp)
# target_link_libraries(darwin PRIVATE covariant_core FFTW3::fftw3f Matplot++::matplot nlohmann_json::nlohmann_json cxxopts::cxxopts Eigen3::Eigen ${MATH_LIB})

add_executable(leonard leonard.cpp Leonard.hpp FlowJo.cpp FlowJo.hpp Pursuit.cpp Pursuit.hpp Reports.cpp Reports.hpp)
target_link_libraries(leonard PRIVATE covariant_core ftxui::screen ftxui::dom ftxui::component LibXml2::LibXml2 LibXslt::LibXslt FFTW3::fftw3f Eigen3::Eigen)

# Standalone gating test executable (using the main() function from Gating.cpp)
# add_executable(gating_test Gating.cpp)
# target_compile_definitions(gating_test PRIVATE BUILD_GATING_TEST)
# target_link_libraries(gating_test PRIVATE covariant_core)

# Set output names to match your Makefile targets
# set_target_properties(riemann PROPERTIES OUTPUT_NAME "riemann")
# set_target_properties(darwin PROPERTIES OUTPUT_NAME "darwin")
set_target_properties(leonard PROPERTIES OUTPUT_NAME "leonard")