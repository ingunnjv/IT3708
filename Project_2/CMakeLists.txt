cmake_minimum_required(VERSION 3.9)
project(Project_2)



set(CMAKE_CXX_STANDARD 11)
INCLUDE_DIRECTORIES(C:/Users/LasseBot/Desktop/eigen-eigen-5a0156e40feb)
#INCLUDE_DIRECTORIES(C:/Program Files/eigen-eigen-5a0156e40feb)
add_executable(Project_2 src/main.cpp src/image_loader.cpp src/image_loader.h)
