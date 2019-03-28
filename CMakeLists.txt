cmake_minimum_required (VERSION 3.1)

project ("RayTracing")

set(ENABLE_UNIT_TESTS CACHE BOOL OFF)
set(CMAKE_CXX_STANDARD 14)

add_subdirectory("src")

if(ENABLE_UNIT_TESTS)
add_subdirectory("tests")
endif(ENABLE_UNIT_TESTS)

include_directories ("data")
