cmake_minimum_required(VERSION 3.15)

project(mpcd)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(PYBIND11_FINDPYTHON ON)

find_package(Python COMPONENTS Interpreter Development REQUIRED)
find_package(pybind11 REQUIRED)
find_package(Eigen3 REQUIRED)

add_library(mpcd_core
    cpp/src/System.cpp
    cpp/src/MPCD.cpp
    cpp/src/Sampling.cpp
    cpp/src/Observables.cpp
)

target_include_directories(mpcd_core PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/cpp/include
)

target_link_libraries(mpcd_core PUBLIC
    Eigen3::Eigen
)

pybind11_add_module(mpcd_cpp
    cpp/bindings/python_bindings.cpp
)

target_include_directories(mpcd_cpp PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/cpp/include
)

target_link_libraries(mpcd_cpp PRIVATE
    mpcd_core
    Eigen3::Eigen
)