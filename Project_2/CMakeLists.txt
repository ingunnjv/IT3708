cmake_minimum_required(VERSION 3.9)
project(Project_2)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -w")

find_package(OpenCV REQUIRED PATHS "C:/OpenCV/opencv/mingw-build")

INCLUDE_DIRECTORIES(C:/Users/Public/eigen-eigen-5a0156e40feb)
#INCLUDE_DIRECTORIES(C:/Program Files/eigen-eigen-5a0156e40feb)

include_directories(${OpenCV_INCLUDE_DIRS})

add_executable(Project_2 src/main.cpp src/image_loader.cpp src/image_loader.h src/genotype.cpp src/genotype.h src/nsga2.cpp src/nsga2.h src/utils.h)

# add libs you need
set(OpenCV_LIBS opencv_core opencv_imgproc opencv_highgui opencv_imgcodecs)
target_link_libraries(Project_2 ${OpenCV_LIBS})

