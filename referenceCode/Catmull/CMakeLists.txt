cmake_minimum_required(VERSION 3.13)
project(catmull)
add_subdirectory(vcglib)
set(CMAKE_CXX_STANDARD 14)
set(SOURCES
	catmull.cpp 
	)

add_executable(catmull
	${SOURCES})

target_link_libraries(
	catmull
	PUBLIC
		vcglib
	)
