cmake_minimum_required(VERSION 3.5)
project(udpclust)

set(CMAKE_CXX_STANDARD 14)
include_directories( src)
set(SOURCE_FILES main.cpp src/cluster.cpp)
add_executable(udpclust ${SOURCE_FILES})