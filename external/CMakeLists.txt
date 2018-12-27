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

# Google test
add_subdirectory(googletest)
set(LIBS ${LIBS} gtest_main)


# Propagate libraries to the root
set(LIBS ${LIBS} PARENT_SCOPE)
