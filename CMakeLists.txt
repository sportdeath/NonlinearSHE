cmake_minimum_required (VERSION 3.4)
project (YASHE)

set(CMAKE_CXX_STANDARD 11)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
  #set(CMAKE_CXX_FLAGS "-Wall -Wextra")
  set(CMAKE_CXX_FLAGS_DEBUG "-g")
  set(CMAKE_CXX_FLAGS_RELEASE "-O3")
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
  ${PROJECT_SOURCE_DIR}/cipherText.cpp
  ${PROJECT_SOURCE_DIR}/numberTheory.cpp
  ${PROJECT_SOURCE_DIR}/functions.cpp
)

include_directories("${PROJECT_INCLUDE_DIR}")

#######################
## External Libraries
#######################

add_subdirectory(external)

####################
## Library Creation
####################

add_library(${PROJECT_NAME} SHARED ${PROJECT_SOURCES})
target_link_libraries(${PROJECT_NAME} ${LIBS})
set(LIBS ${LIBS} ${PROJECT_NAME})

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
MAKE_TEST(timing_test)

######################
## Building examples
######################

option(BUILD_EXAMPLES "BUILD_EXAMPLES" OFF)
if(BUILD_EXAMPLES)
  add_subdirectory(examples)
endif(BUILD_EXAMPLES)
