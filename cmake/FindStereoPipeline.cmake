# Find Stereo Pipeline
#
# == Using Header-Only libraries from within Stereo Pipeline: ==
#
#  find_package( StereoPipeline 2.0 )
#  if(StereoPipeline_FOUND)
#     include_directories(${StereoPipeline_INCLUDE_DIRS})
#     add_executable(foo foo.cc)
#  endif()
#
# == Using actual libraries from within StereoPipeline: ==
#
#  find_package( StereoPipeline 2.0 COMPONENTS core isisio ... )
#
#  if(StereoPipeline_FOUND)
#     include_directories(${StereoPipeline_INCLUDE_DIRS})
#     add_executable(foo foo.cc)
#     target_link_libraries(foo ${StereoPipeline_LIBRARIES})
#  endif()
#
# ===========================================================
#
#  Copyright (c) 2010 Zachary Moratto
#
#  Redistribution AND use is allowed according to the terms
#  of the NEW BSD license.
#

# This macro must define:
#  StereoPipeline_FOUND             < Conditional
#  StereoPipeline_INCLUDE_DIR       < Paths
#  StereoPipeline_LIBRARIES         < Paths as well

# Find Path



find_path(StereoPipeline_DIR
  NAMES  include/asp/Core.h
         lib/libaspCore.dylib
  HINTS  ${StereoPipeline_INCLUDE_DIR}
)

find_path(StereoPipeline_INCLUDE_DIR
  NAMES  asp/Core.h
  HINTS  ${StereoPipeline_DIR}
  PATH_SUFFIXES include
)

# Searching for each library that makes up a component
foreach(COMPONENT ${StereoPipeline_FIND_COMPONENTS})
  string(TOUPPER ${COMPONENT} UPPERCOMPONENT)
  set( StereoPipeline_${UPPERCOMPONENT}_LIBRARY "StereoPipeline_${UPPERCOMPONENT}_LIBRARY-NOTFOUND")

  string(REGEX REPLACE "^(.)(.+)$" "\\1" first_half "${UPPERCOMPONENT}")
  string(REGEX REPLACE "^(.)(.+)$" "\\2" second_half "${COMPONENT}")

  set(LIBRARY_NAME "${first_half}${second_half}")
  if ( ${UPPERCOMPONENT} STREQUAL "ISISIO" )
    set(LIBRARY_NAME "IsisIO")
  elseif( ${UPPERCOMPONENT} STREQUAL "PHOTOMETRYTK" )
    set(LIBRARY_NAME "PhotometryTK")
  endif()

  find_library( StereoPipeline_${UPPERCOMPONENT}_LIBRARY
    NAMES  asp${LIBRARY_NAME}
    HINTS  ${StereoPipeline_DIR}
    PATH_SUFFIXES lib lib64
    )

  if(StereoPipeline_${UPPERCOMPONENT}_LIBRARY)
    set(StereoPipeline_${UPPERCOMPONENT}_FOUND TRUE CACHE INTERNAL "If the ASP ${UPPERCOMPONENT} library was found")
  endif()

endforeach(COMPONENT)

# Deciding if ASP was found
set(StereoPipeline_INCLUDE_DIRS ${StereoPipeline_INCLUDE_DIR})

if(StereoPipeline_INCLUDE_DIR)
  set(StereoPipeline_FOUND TRUE)
else(StereoPipeline_INCLUDE_DIR)
  set(StereoPipeline_FOUND FALSE)
endif(StereoPipeline_INCLUDE_DIR)

# Closing Messages
if(StereoPipeline_FOUND)
  message(STATUS "Found the following StereoPipeline libraries:")
  foreach( COMPONENT ${StereoPipeline_FIND_COMPONENTS})
    string(TOUPPER ${COMPONENT} UPPERCOMPONENT )
    if ( StereoPipeline_${UPPERCOMPONENT}_FOUND )
      message(STATUS "  ${COMPONENT}")
      set(StereoPipeline_LIBRARIES ${StereoPipeline_LIBRARIES} ${StereoPipeline_${UPPERCOMPONENT}_LIBRARY})
    endif()
  endforeach()
else(StereoPipeline_FOUND)
  message(WARNING "Unable to find requested StereoPipeline libraries")
endif(StereoPipeline_FOUND)


