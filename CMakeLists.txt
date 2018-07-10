cmake_minimum_required(VERSION 2.8.3)
project(eclib)

set(CMAKE_CXX_FLAGS "-std=c++14 -O3 ${CMAKE_CXX_FLAGS}")

INCLUDE_DIRECTORIES( 
	include 
)

add_library(eclib SHARED 
	src/CliffordLib.cpp
)

