cmake_minimum_required(VERSION 3.0.0...3.5)
project(polytopewalk)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "-O2 -fPIC")

find_package(Python COMPONENTS Interpreter Development.Module REQUIRED)
find_package(pybind11 CONFIG REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(ifopt REQUIRED)

set(POLYTOPEWALK_DIR ${CMAKE_CURRENT_SOURCE_DIR})
set(SOURCES_DIR ${POLYTOPEWALK_DIR}/src)
set(TESTS_DIR ${POLYTOPEWALK_DIR}/tests)

add_subdirectory("src")
# set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})
#add_subdirectory ("docs")


pybind11_add_module(${PROJECT_NAME} MODULE pybind11/PybindExt.cpp)
target_link_libraries(${PROJECT_NAME} PUBLIC randomwalk pybind11::module Eigen3::Eigen ifopt::ifopt_ipopt)
target_include_directories(${PROJECT_NAME} PRIVATE ${Eigen3_INCLUDE_DIRS})
install(TARGETS ${PROJECT_NAME} LIBRARY DESTINATION .)

pybind11_add_module(sparsepolytopewalk MODULE pybind11/SparsePybindExt.cpp)
target_link_libraries(sparsepolytopewalk PUBLIC sparserandomwalk pybind11::module Eigen3::Eigen ifopt::ifopt_ipopt)
target_include_directories(sparsepolytopewalk PRIVATE ${Eigen3_INCLUDE_DIRS})
install(TARGETS sparsepolytopewalk LIBRARY DESTINATION .)

add_subdirectory("tests")