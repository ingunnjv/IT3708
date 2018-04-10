cmake_minimum_required(VERSION 3.9)
project(Project_2)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -w")

INCLUDE_DIRECTORIES(C:/Users/Public/eigen-eigen-5a0156e40feb)

add_executable(Project_3 src/main.cpp src/jssp.cpp src/jssp.h)



