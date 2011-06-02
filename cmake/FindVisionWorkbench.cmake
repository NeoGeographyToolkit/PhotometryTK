# Find Vision Workbench
#
#  Copyright (c) 2009-2011 Zachary Moratto
#  Large parts of this work we're inspired by FindBoost
#
#  Redistribution AND use is allowed according to the terms
#  of the NEW BSD license.
#
# == Using Header-Only libraries from within Vision Workbench: ==
#
#  find_package( VisionWorkbench )
#  if(VISIONWORKBENCH_FOUND)
#     include_directories(${VISIONWORKBENCH_INCLUDE_DIRS})
#     add_executable(foo foo.cc)
#  endif()
#
# == Using actual libraries from within Vision Workbench: ==
#
#  find_package( VisionWorkbench COMPONENTS core image fileio camera ... )
#
#  if(VISIONWORKBENCH_FOUND)
#     include_directories(${VISIONWORKBENCH_INCLUDE_DIRS})
#     add_executable(foo foo.cc)
#     target_link_libraries(foo ${VISIONWORKBENCH_LIBRARIES})
#  endif()
#
# ===========================================================
#
# Variables used by this module, they can change the default behaviour and
# need to be set before calling find_package:
#
#  VISIONWORKBENCH_USE_STATIC_LIBS    Can be set to ON to force the use of the
#                                     static libraries. Default is OFF.
#
#  VISIONWORKBENCH_ROOT               The preferred installation prefix for searching
#                                     for Vision Workbench.
#
# ===========================================================
#
# Variables defined by this module:
#
#  VISIONWORKBENCH_FOUND              System has Vision Workbench with options requested
#
#  VISIONWORKBENCH_INCLUDE_DIRS       Vision Workbench include directories
#
#  VISIONWORKBENCH_LIBRARIES          Link to these to use the Vision Workbench libraries
#
#  VISIONWORKBENCH_LIBRARY_DIRS       The path to where the Vison Workbench library files are.
#
# For each component you specify in find_package(), the folloing (UPPER-CASE)
# variables are set. You can use these variables if you would like to pick and
# choose components for your targets instead of just using VISIONWORKBENCH_LIBRARIES
#
#  VISIONWORKBENCH_${COMPONENT}_FOUND    True if the Vision Workbench "component" was found
#
#  VISIONWORKBENCH_${COMPONENT}_LIBRARY  Contains the libraries for the specified
#                                        VisionWorkbench "component".

find_path(VISIONWORKBENCH_INCLUDE_DIR
  NAMES  vw/vw.h
  HINTS  $ENV{VISIONWORKBENCH_ROOT}
         ${VISIONWORKBENCH_ROOT}
         $ENV{HOME}/projects/VisionWorkbench
  PATH_SUFFIXES include
)

# Check if we are only looking for static
if(VISIONWORKBENCH_USE_STATIC_LIBS)
  if(WIN32)
    set(CMAKE_FIND_LIBRARY_SUFFIXES .lib .a ${CMAKE_FIND_LIBRARY_SUFFIXES})
  else(WIN32)
    set(CMAKE_FIND_LIBRARY_SUFFIXES .a)
  endif(WIN32)
endif(VISIONWORKBENCH_USE_STATIC_LIBS)

# Searching for each library that makes up a component
foreach(COMPONENT ${VisionWorkbench_FIND_COMPONENTS})
  string(TOUPPER ${COMPONENT} UPPERCOMPONENT)
  set( VISIONWORKBENCH_${UPPERCOMPONENT}_LIBRARY "VISIONWORKBENCH_${UPPERCOMPONENT}_LIBRARY-NOTFOUND"
    CACHE FILEPATH "The VisionWorkbench ${COMPONENT} library")

  string(REGEX REPLACE "^(.)(.+)$" "\\1" first_half "${UPPERCOMPONENT}")
  string(REGEX REPLACE "^(.)(.+)$" "\\2" second_half "${COMPONENT}")

  set(LIBRARY_NAME "${first_half}${second_half}")
  if ( ${UPPERCOMPONENT} STREQUAL "VW" )
    set(LIBRARY_NAME "")
  elseif( ${UPPERCOMPONENT} STREQUAL "INTERESTPOINT" )
    set(LIBRARY_NAME "InterestPoint")
  elseif( ${UPPERCOMPONENT} STREQUAL "HDR" )
    set(LIBRARY_NAME "HDR")
  elseif( ${UPPERCOMPONENT} STREQUAL "FILEIO")
    set(LIBRARY_NAME "FileIO")
  elseif( ${UPPERCOMPONENT} STREQUAL "BUNDLEADJUSTMENT")
    set(LIBRARY_NAME "BundleAdjustment")
  endif()

  find_library( VISIONWORKBENCH_${UPPERCOMPONENT}_LIBRARY
    NAMES  vw${LIBRARY_NAME}
    HINTS  ${VISIONWORKBENCH_INCLUDE_DIR}/..
           ${VISIONWORKBENCH_ROOT}
           $ENV{VISIONWORKBENCH_ROOT}
           $ENV{HOME}/projects/VisionWorkbench
    PATH_SUFFIXES lib lib64
    )
  mark_as_advanced( VISIONWORKBENCH_${UPPERCOMPONENT}_LIBRARY )

  if(VISIONWORKBENCH_${UPPERCOMPONENT}_LIBRARY)
    set(VISIONWORKBENCH_${UPPERCOMPONENT}_FOUND TRUE CACHE INTERNAL "If the VW ${UPPERCOMPONENT} library was found")
  endif()

endforeach(COMPONENT)

# Deciding if VW was found
set(VISIONWORKBENCH_INCLUDE_DIRS ${VISIONWORKBENCH_INCLUDE_DIR})

if(VISIONWORKBENCH_INCLUDE_DIR)
  set(VISIONWORKBENCH_FOUND TRUE)
else(VISIONWORKBENCH_INCLUDE_DIR)
  set(VISIONWORKBENCH_FOUND FALSE)
endif(VISIONWORKBENCH_INCLUDE_DIR)

# Closing Messages
if(VISIONWORKBENCH_FOUND)
  message(STATUS "Found the following VisionWorkbench libraries:")
  foreach( COMPONENT ${VisionWorkbench_FIND_COMPONENTS})
    string(TOUPPER ${COMPONENT} UPPERCOMPONENT )
    if ( VISIONWORKBENCH_${UPPERCOMPONENT}_FOUND )
      message(STATUS "  ${COMPONENT}")
      set(VISIONWORKBENCH_LIBRARIES ${VISIONWORKBENCH_LIBRARIES} ${VISIONWORKBENCH_${UPPERCOMPONENT}_LIBRARY})
    endif()
  endforeach()
else(VISIONWORKBENCH_FOUND)
  message(SEND_ERROR "Unable to find requested VisionWorkbench libraries")
endif(VISIONWORKBENCH_FOUND)