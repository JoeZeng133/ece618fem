project(ece618)
cmake_minimum_required (VERSION 3.9)
set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_STANDARD 11)

#However, the file(GLOB...) allows for wildcard additions:
file(GLOB SOURCES "src/*.cpp")

add_executable(ece618 ${SOURCES})