project (FHEImageProcessing)

set(PROJECT_INCLUDE_DIR "${PROJECT_SOURCE_DIR}/include")
set(PROJECT_SOURCE_DIR "${PROJECT_SOURCE_DIR}/lib")
set(PROJECT_SOURCES
  ${PROJECT_SOURCE_DIR}/imageFunctions.cpp
)
include_directories("${PROJECT_INCLUDE_DIR}")

# CImg
find_package(CIMG REQUIRED)
include_directories(${CIMG_INCLUDE_DIRS})
link_directories(${cimg_system_libs_dir})
set(LIBS ${LIBS} ${CIMG_SYSTEM_LIBS})
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CIMG_CFLAGS}")
if(NOT CIMG_FOUND)
  message(ERROR " CIMG not found!")
endif(NOT CIMG_FOUND)

add_executable(imageTransform
               src/imageTransform.cpp
               ${PROJECT_SOURCES}
)
target_link_libraries(
  imageTransform
  ${LIBS}
)

add_executable(RGBtoYCbCr
               src/RGBtoYCbCr.cpp
               ${PROJECT_SOURCES}
)
target_link_libraries(
  RGBtoYCbCr
  ${LIBS}
)

file(COPY resources DESTINATION "${CMAKE_CURRENT_BINARY_DIR}")
