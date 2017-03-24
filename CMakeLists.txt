cmake_minimum_required(VERSION 3.5)

project(xfem)

file(GLOB FILES *.cpp)

ADD_DEFINITIONS(-openmp)

include_directories(include/gmsh)
include_directories(include/gmm)

link_directories(lib)
set(EXECUTABLE_OUTPUT_PATH ../bin)


add_executable(xfem ${FILES})

target_link_libraries(xfem Gmsh lapack blas)
