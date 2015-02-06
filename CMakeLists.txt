
cmake_minimum_required( VERSION 2.6 )
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

if (NOT CMAKE_BUILD_TYPE)
	message(STATUS "No build type selected, default to Release")
	set(CMAKE_BUILD_TYPE "Release")
endif()


project( Polynomial )

find_package( Eigen3 REQUIRED )
include_directories( ${EIGEN3_INCLUDE_DIR} )

include_directories( /usr/local/include )

include_directories( $(Polynomial) )

add_executable( testpolynomial testpolynomial.cpp )

