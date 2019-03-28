cmake_minimum_required (VERSION 3.1)

set(APPLICATION_NAME Raytracing)

message(STATUS "Configuring ${APPLICATION_NAME}")

find_package( Boost 1.65 COMPONENTS program_options filesystem REQUIRED )
include_directories( ${Boost_INCLUDE_DIR} )
add_executable (${APPLICATION_NAME} "Raytracing.cpp" "Utilities.h" "Utilities.cpp")
target_link_libraries(${APPLICATION_NAME} ${Boost_LIBRARIES} )
