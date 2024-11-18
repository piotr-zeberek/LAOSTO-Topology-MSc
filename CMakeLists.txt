cmake_minimum_required(VERSION 3.10.0)
project (system2D)

# prevent for in-source building
set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)
set(CMAKE_DISABLE_SOURCE_CHANGES  ON)
if ("${CMAKE_SOURCE_DIR}" STREQUAL "${CMAKE_BINARY_DIR}")
  message(SEND_ERROR "In-source builds are not allowed.")
endif ()

add_definitions(-DPROJECTPATH=\"${PROJECT_SOURCE_DIR}\")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")

set(INCLUDES ${PROJECT_SOURCE_DIR}/include)
include_directories(${INCLUDES})

file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/data)
file(GLOB SOURCES "${PROJECT_SOURCE_DIR}/src/*.cpp")

add_executable(${PROJECT_NAME} main.cpp ${SOURCES})

# include Eigen
find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

# Link OpenMP
find_package(OpenMP REQUIRED)
target_link_libraries(${PROJECT_NAME} PUBLIC OpenMP::OpenMP_CXX)

# Link Armadillo
find_package(Armadillo REQUIRED)
include_directories(${ARMADILLO_INCLUDE_DIRS})
target_link_libraries(${PROJECT_NAME} PUBLIC ${ARMADILLO_LIBRARIES})

# Find the ARPACK library
find_library(ARPACK_LIB NAMES arpack PATHS /usr/lib /usr/local/lib)

# Add the ARPACK library to the target link libraries
target_link_libraries(${PROJECT_NAME} PUBLIC ${ARPACK_LIB})

# Find the SuperLU library
find_library(SUPERLU_LIB NAMES superlu PATHS /usr/lib /usr/local/lib)

# Add the SuperLU library to the target link libraries
target_link_libraries(${PROJECT_NAME} PUBLIC ${SUPERLU_LIB})

# speed up
target_compile_options(${PROJECT_NAME} PUBLIC -march=native)
