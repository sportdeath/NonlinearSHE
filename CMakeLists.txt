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

include(ExternalProject)

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

# Google Test
# For testing
ExternalProject_Add(
  gtest
  GIT_REPOSITORY https://github.com/google/googletest.git
  PREFIX ${CMAKE_CURRENT_BINARY_DIR}/gtest
  TIMEOUT 10
  UPDATE_COMMAND ""
  INSTALL_COMMAND "")

ExternalProject_Get_Property(gtest source_dir binary_dir)
include_directories("${source_dir}/googletest/include")
set(LIBS ${LIBS} "${binary_dir}/googlemock/gtest/libgtest.a")

##################
## Unit tests
#################

add_executable(functions_test 
       unitTests/functions_test.cpp
       source/functions.cpp
       source/functions.hpp
       )

add_dependencies(functions_test gtest)

target_link_libraries(
  functions_test
  ${LIBS}
)

add_executable(numberTheory_test
         unitTests/numberTheory_test.cpp
         source/numberTheory.cpp
         source/numberTheory.hpp
)

add_dependencies(numberTheory_test gtest)

target_link_libraries(
  numberTheory_test
  ${LIBS}
)

add_executable(yashe_test
         unitTests/yashe_test.cpp
         source/numberTheory.hpp
         source/numberTheory.cpp
         source/yashe.cpp
         source/yashe.hpp
         source/cipherText.cpp
         source/cipherText.hpp
         source/functions.cpp
         source/functions.hpp
         serialization/ZZ.hpp
         serialization/ZZ_p.hpp
         serialization/ZZX.hpp
         serialization/ZZ_pX.hpp
)

add_dependencies(yashe_test gtest)

target_link_libraries(
  yashe_test
  ${LIBS}
)
