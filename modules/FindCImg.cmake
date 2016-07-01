# - Try to find CImg lib
#
# The following variables are defined
#
#  CIMG_FOUND - system has CImg lib
#  CIMG_INCLUDE_DIRS - the CImg include directory
#  CIMG_SYSTEM_LIBS - external libraries that CImg uses
#  CIMG_SYSTEM_LIBS_DIR - external library directories
#  CIMG_CFLAGS - compilation flags


if (NOT CIMG_INCLUDE_DIR)
  find_path(CIMG_INCLUDE_DIR
    NAMES CImg.h
    PATHS
      ${CMAKE_INSTALL_PREFIX}/include
      /usr/include
  )
  mark_as_advanced(CIMG_INCLUDE_DIR)
endif(NOT CIMG_INCLUDE_DIR)
if (NOT CIMG_INCLUDE_DIR-NOTFOUND)
  set(CIMG_FOUND true)
endif(NOT CIMG_INCLUDE_DIR-NOTFOUND)
list(APPEND CIMG_INCLUDE_DIRS
  ${CIMG_INCLUDE_DIR}
)

# ### CIMG related stuff
# Flags to enable fast image display, using the XSHM library.
SET(CIMG_XSHM_CCFLAGS  -Dcimg_use_xshm)

# Flags to enable screen mode switching, using the XRandr library.
SET(CIMG_XRANDR_CCFLAGS  -Dcimg_use_xrandr)

# Flags to enable native support for JPEG image files, using the JPEG library.
# ( http://www.ijg.org/ )
SET(CIMG_JPEG_CCFLAGS  -Dcimg_use_jpeg)

# Flags to enable native support for TIFF image files, using the TIFF library.
# ( http://www.libtiff.org/ )
SET(CIMG_TIFF_CCFLAGS  -Dcimg_use_tiff)

# Flags to enable native support for PNG image files, using the PNG library.
# ( http://www.libpng.org/ )
SET(CIMG_PNG_CCFLAGS  -Dcimg_use_png)

# ### Search Additional Libraries ##########
FIND_PACKAGE(JPEG)
FIND_PACKAGE(TIFF)
FIND_PACKAGE(PNG)

if(NOT WIN32)
  FIND_PACKAGE(X11)
  FIND_PACKAGE(Threads REQUIRED)
endif()

# #### End of additional libraries search ##########

### Configure Paths according to detected packages
if(TIFF_FOUND)
  get_filename_component(TIFF_LIB_DIRS ${TIFF_LIBRARIES} PATH)
  SET(CIMG_CFLAGS "${CIMG_CFLAGS} ${CIMG_TIFF_CCFLAGS}")
  link_directories(${TIFF_LIB_DIRS})
  include_directories(${TIFF_INCLUDE_DIR})
  SET(CIMG_SYSTEM_LIBS ${CIMG_SYSTEM_LIBS} ${TIFF_LIBRARIES})
endif()

if(JPEG_FOUND)
  get_filename_component(JPEG_LIB_DIRS ${JPEG_LIBRARIES} PATH)
  SET(CIMG_CFLAGS "${CIMG_CFLAGS} ${CIMG_JPEG_CCFLAGS}")
  link_directories(${JPEG_LIB_DIRS})
  include_directories(${JPEG_INCLUDE_DIR})
  SET(CIMG_SYSTEM_LIBS ${CIMG_SYSTEM_LIBS} ${JPEG_LIBRARIES})
endif()

if(NOT APPLE)
  if(NOT WIN32)
    if(X11_FOUND)
      SET(CIMG_CFLAGS "${CIMG_CFLAGS} ${CIMG_XSHM_CCFLAGS} ${CIMG_XRANDR_CCFLAGS}")
      SET(CIMG_SYSTEM_LIBS ${CIMG_SYSTEM_LIBS} Xext Xrandr)
    endif()
  endif(NOT WIN32)
endif(NOT APPLE)

if(X11_FOUND)
  link_directories(${X11_LIB_DIRS})
  include_directories(${X11_INCLUDE_DIR})
  SET( CIMG_SYSTEM_LIBS ${CIMG_SYSTEM_LIBS} ${X11_LIBRARIES} )
endif()

if (NOT WIN32)
  SET( CIMG_SYSTEM_LIBS ${CIMG_SYSTEM_LIBS} ${CMAKE_THREAD_LIBS_INIT} )
endif()

if( WIN32)
  SET( CIMG_SYSTEM_LIBS  ${CIMG_SYSTEM_LIBS}  gdi32 )
endif()

# Add CIMG Flags to Compilation Flags
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CIMG_CFLAGS}")

foreach(program ${CIMG_FILES})
  add_executable(${program} ${program}.cpp)
  target_link_libraries(${program} ${CIMG_SYSTEM_LIBS} )
endforeach(program)
