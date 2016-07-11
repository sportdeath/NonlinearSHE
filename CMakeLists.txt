cmake_minimum_required (VERSION 3.4)
project (YASHE)

include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
if(COMPILER_SUPPORTS_CXX11)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
else()
  message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no c++11 support.")
endif()

#Adding cmake modules
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/modules/)

#######################
## Internal Libraries
#######################

set(PROJECT_INCLUDE_DIR "${PROJECT_SOURCE_DIR}/include")
set(PROJECT_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src")
set(PROJECT_SOURCES
  ${PROJECT_SOURCE_DIR}/yashe.cpp
  ${PROJECT_SOURCE_DIR}/ciphertext.cpp
  ${PROJECT_SOURCE_DIR}/numberTheory.cpp
  ${PROJECT_SOURCE_DIR}/functions.cpp
)

include_directories("${PROJECT_INCLUDE_DIR}")

#######################
## External Libraries
#######################

# NTL
find_package(NTL REQUIRED)
include_directories(${NTL_INCLUDE_DIR})
set(LIBS ${LIBS} ${NTL_LIBRARY})
if(NOT NTL_FOUND)
  message(ERROR " NTL not found!")
endif(NOT NTL_FOUND)

# GMP
find_package(GMP REQUIRED)
include_directories(${GMP_INCLUDE_DIR})
set(LIBS ${LIBS} ${GMP_LIBRARIES})
if(NOT GMP_FOUND)
  message(ERROR " GMP not found!")
endif(NOT GMP_FOUND)

# Boost
# For serialization (writing class to file)
find_package(Boost COMPONENTS serialization REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})
link_directories(${Boost_LIBRARY_DIRS})
set(LIBS ${LIBS} ${Boost_LIBRARIES})
if(NOT Boost_FOUND)
  message(ERROR " Boost not found!")
endif(NOT Boost_FOUND)

####################
## Library Creation
####################

add_library(${PROJECT_NAME} SHARED ${PROJECT_SOURCES})
target_link_libraries(${PROJECT_NAME} ${LIBS})
set(LIBS ${LIBS} ${PROJECT_NAME})

##################
## Google Test
#################

# Google Test
# For testing
include(ExternalProject)
ExternalProject_Add(
  gtest
  GIT_REPOSITORY https://github.com/google/googletest.git
  PREFIX ${CMAKE_CURRENT_BINARY_DIR}/gtest
  TIMEOUT 10
  UPDATE_COMMAND ""
  INSTALL_COMMAND "")

enable_testing()

ExternalProject_Get_Property(gtest source_dir binary_dir)
include_directories("${source_dir}/googletest/include")
set(LIBS ${LIBS} "${binary_dir}/googlemock/gtest/libgtest.a")

#####################################
## Unit Tests with CTest and GTest
#####################################

macro(MAKE_TEST test_name)
  add_executable(${test_name}
         test/${test_name}.cpp
         )
  target_link_libraries(${test_name} ${LIBS})
  add_test(${test_name} ${test_name})
endmacro()

MAKE_TEST(functions_test)
MAKE_TEST(numberTheory_test)
MAKE_TEST(yashe_test)

######################
## Building examples
######################

option(BUILD_EXAMPLES "BUILD_EXAMPLES" OFF)
if(BUILD_EXAMPLES)
  add_subdirectory(examples)
endif(BUILD_EXAMPLES)
