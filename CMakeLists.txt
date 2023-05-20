cmake_minimum_required(VERSION 3.13)
project(polytopewalk)

set(CMAKE_OSX_ARCHITECTURES "arm64")
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "-O2")

find_package(Python COMPONENTS Interpreter Development.Module REQUIRED)

find_package(ifopt REQUIRED)
find_package(pybind11 CONFIG REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(Catch2 REQUIRED)

# Include the local folder "include" for header files
include_directories(include)

# Define the source files for the library
file(GLOB SOURCES "src/*.cpp")
list(FILTER SOURCES EXCLUDE REGEX ".*PybindExt\\.cpp$")
# message("The value of SOURCES is: ${SOURCES}")

pybind11_add_module(${PROJECT_NAME} MODULE src/PybindExt.cpp ${SOURCES})
target_include_directories(${PROJECT_NAME} PRIVATE ${Eigen3_INCLUDE_DIRS})
target_link_libraries(${PROJECT_NAME} PUBLIC pybind11::module Eigen3::Eigen ifopt::ifopt_ipopt)

add_executable(tests tests/test.cpp ${SOURCES})
target_include_directories(tests PRIVATE ${Eigen3_INCLUDE_DIRS})
target_link_libraries(tests PUBLIC pybind11::module Eigen3::Eigen ifopt::ifopt_ipopt Catch2::Catch2WithMain)

install(TARGETS ${PROJECT_NAME}
        LIBRARY DESTINATION .
)