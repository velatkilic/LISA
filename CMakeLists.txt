cmake_minimum_required(VERSION 3.19)
project(LISA)

add_subdirectory(extern/pybind11)
include_directories(extern/Eigen)

option(BUILD_PIP "True when pip installing" False)

# lisa is a pybind11 module if building pip
# if debugging with visual studio set BUILD_PIP to False
if(BUILD_PIP)
	pybind11_add_module(pylisa "")
	add_subdirectory(wrappers)
else()
	add_library(pylisa "")
endif()

add_subdirectory("include")
add_subdirectory(src)

# Specify include directories
target_include_directories(pylisa PUBLIC
	"$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
	"$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>"
)

# link multi-threading library
find_package(Threads)
target_link_libraries(pylisa PRIVATE ${CMAKE_THREAD_LIBS_INIT})

# cpp tests for when pip is not built
if(NOT BUILD_PIP)
	add_subdirectory(tests)
endif()