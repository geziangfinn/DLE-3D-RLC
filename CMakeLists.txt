cmake_minimum_required(VERSION 3.14)
project(ZST_DME)

set(CMAKE_CXX_STANDARD 17)
#add_definitions(-std=c99)
#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c90 -pedantic")


add_executable(ZST_DME dle3dmain.cpp DME.h OutputFiles.h DLE.h)
