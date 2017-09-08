cmake_minimum_required(VERSION 3.5)

project(xfem CXX C Fortran)

file(GLOB FILES *.cpp *.f)

find_package(LAPACK REQUIRED)
add_definitions(-DGMM_USES_LAPACK)

include_directories(include/gmsh)
include_directories(include/gmm)

link_directories(lib)
set(EXECUTABLE_OUTPUT_PATH ../bin)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

add_executable(xfem ${FILES})

target_link_libraries(xfem Gmsh ${MUMPS_LIBRARIES} ${LAPACK_LIBRARIES} )
