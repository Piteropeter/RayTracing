cmake_minimum_required (VERSION 3.1)

set(APPLICATION_NAME raytracing)

message(STATUS "Configuring ${APPLICATION_NAME}")

FIND_PACKAGE( Boost 1.68 COMPONENTS program_options filesystem REQUIRED )
#add_library(Reader "Reader.cpp" "Reader.h")
#add_library(Solver "Solver.cpp" "Solver.h" "Utilities.h" "BruteForce.cpp" "BruteForce.h" "BranchAndBound.cpp" "BranchAndBound.h" "Dynamic.cpp" "Dynamic.h" "Tabu.cpp" "Tabu.h" "Genetic.cpp" "Genetic.h")
INCLUDE_DIRECTORIES( ${Boost_INCLUDE_DIR} )
add_executable (${APPLICATION_NAME} "raytracing.cpp" "Utilities.h" "Utilities.cpp")
TARGET_LINK_LIBRARIES(${APPLICATION_NAME} ${Boost_LIBRARIES} )

target_link_libraries(${APPLICATION_NAME})
