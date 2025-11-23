cmake_minimum_required (VERSION 3.5)

project (blackray)

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -O3 -fopenmp") #-Werror
set (src_dir_rt "${PROJECT_SOURCE_DIR}/src/raytracer/")
set (src_dir_conv "${PROJECT_SOURCE_DIR}/src/convolver/")

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/")

file (GLOB source_files_conv "${src_dir_conv}/convolver.cpp")
file (GLOB source_files_rt "${src_dir_rt}/main.cpp")

add_executable (rt ${source_files_rt})
add_executable (conv ${source_files_conv})
