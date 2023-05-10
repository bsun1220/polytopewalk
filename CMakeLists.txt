cmake_minimum_required(VERSION 3.13)
project(polytopewalk)
set(CMAKE_OSX_ARCHITECTURES "x86_64")
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "-O2")

find_package(Python COMPONENTS Interpreter Development.Module REQUIRED)

find_package(ifopt REQUIRED)
find_package(pybind11 CONFIG REQUIRED)
find_package(Eigen3 REQUIRED)

# Include the local folder "include" for header files
include_directories(include)

#set(PYTHON_EXECUTABLE /Users/yc443/opt/anaconda3/bin/python)

# Define the source files for the library
file(GLOB SOURCES "src/*.cpp")
#list(FILTER SOURCES EXCLUDE REGEX ".*main\\.cpp$")
#list(FILTER SOURCES EXCLUDE REGEX ".*test\\.cpp$")
# message("The value of SOURCES is: ${SOURCES}")

add_executable(${PROJECT_NAME} test.cpp ${SOURCES})
target_link_libraries(${PROJECT_NAME} PUBLIC pybind11::module Eigen3::Eigen ifopt::ifopt_ipopt)
target_include_directories(${PROJECT_NAME} PRIVATE ${Eigen3_INCLUDE_DIRS})
#pybind11_add_module(${PROJECT_NAME} MODULE src/main.cpp ${SOURCES})
#target_link_libraries(${PROJECT_NAME} PUBLIC pybind11::module Eigen3::Eigen ifopt::ifopt_ipopt)
#target_include_directories(${PROJECT_NAME} PRIVATE ${Eigen3_INCLUDE_DIRS})

# add_executable(mytest1 src/test.cpp ${SOURCES})
# target_link_libraries(mytest1 PRIVATE Eigen3::Eigen ifopt::ifopt_ipopt)
# target_include_directories(mytest1 PRIVATE ${Eigen3_INCLUDE_DIRS})

install(TARGETS ${PROJECT_NAME}
        LIBRARY DESTINATION .
)