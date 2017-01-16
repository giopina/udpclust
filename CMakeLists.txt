cmake_minimum_required(VERSION 3.5)
project(udpclust)

set(CMAKE_CXX_STANDARD 14)
ADD_SUBDIRECTORY(libs/pybind11)
set(PYTHON_INCLUDE_DIRS /home/marscher/anaconda/envs/py35/include/python3.5m/)
include_directories(src ${PYBIND11_INCLUDE_DIR} ${PYTHON_INCLUDE_DIRS})

set(SOURCE_FILES src/cluster.cpp
        src/module.cpp
        )
add_executable(udpclust ${SOURCE_FILES})
#TARGET_LINK_LIBRARIES(LINK_PUBLIC ${PYTHON_LIBRARIES} )