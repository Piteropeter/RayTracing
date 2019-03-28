cmake_minimum_required (VERSION 3.1)

set(APPLICATION_NAME ut)

message(STATUS "Configuring ${APPLICATION_NAME}")

enable_testing()
find_package(gtest REQUIRED)
FIND_PACKAGE(Boost 1.65 COMPONENTS program_options filesystem REQUIRED)
include_directories(${GTEST_INCLUDE_DIR})
INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIR})
include_directories(${CMAKE_SOURCE_DIR}/src)
add_executable(${APPLICATION_NAME} "UT.cpp")
#target_link_libraries(Unit_tests Reader)
#target_link_libraries(Unit_tests Solver)
TARGET_LINK_LIBRARIES(${APPLICATION_NAME} ${GTEST_BOTH_LIBRARIES})
TARGET_LINK_LIBRARIES(${APPLICATION_NAME} ${Boost_LIBRARIES})

add_test(${APPLICATION_NAME} UT)