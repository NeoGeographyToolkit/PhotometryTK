
if (MKL_LIBRARIES)
  set(MKL_FIND_QUIETLY TRUE)
endif (MKL_LIBRARIES)

find_path(MKL_INCLUDE_DIR
  NAMES mkl.h
  HINTS $ENV{MKL_ROOT}
        ${MKL_ROOT}
        /opt/intel/mkl
        /usr/intel/mkl
  PATH_SUFFIXES include
  )

if (NOT MKL_ROOT)
  set(MKL_ROOT ${MKL_INCLUDE_DIR}/..)
endif()

if(CMAKE_MINOR_VERSION GREATER 4)

  if(${CMAKE_HOST_SYSTEM_PROCESSOR} STREQUAL "x86_64")

    find_library(MKL_LIBRARIES
      NAMES mkl_core
      PATHS $ENV{MKL_ROOT}/lib/em64t
            ${MKL_ROOT}/lib/em64t
            /opt/intel/mkl/*/lib/em64t
            /opt/intel/Compiler/*/*/mkl/lib/em64t
            ${LIB_INSTALL_DIR}
            )

    find_library(MKL_GUIDE
      NAMES guide
      PATHS $ENV{MKL_ROOT}/lib/em64t
            ${MKL_ROOT}/lib/em64t
            /opt/intel/mkl/*/lib/em64t
            /opt/intel/Compiler/*/*/mkl/lib/em64t
            /opt/intel/Compiler/*/*/lib/intel64
            ${LIB_INSTALL_DIR}
            )

    if(MKL_LIBRARIES AND MKL_GUIDE)
      set(MKL_LIBRARIES ${MKL_LIBRARIES} mkl_intel_lp64 mkl_sequential ${MKL_GUIDE} pthread)
    endif()

  else(${CMAKE_HOST_SYSTEM_PROCESSOR} STREQUAL "x86_64")

    find_library(MKL_LIBRARIES
      NAMES mkl_core
      PATHS $ENV{MKL_ROOT}/lib/32
            ${MKL_ROOT}/lib/32
            /opt/intel/mkl/*/lib/32
            /opt/intel/Compiler/*/*/mkl/lib/32
            ${LIB_INSTALL_DIR}
            )

    find_library(MKL_GUIDE
      NAMES guide
      PATHS $ENV{MKL_ROOT}/lib/32
            ${MKL_ROOT}/lib/32
            /opt/intel/mkl/*/lib/32
            /opt/intel/Compiler/*/*/mkl/lib/32
            /opt/intel/Compiler/*/*/lib/intel32
            ${LIB_INSTALL_DIR}
            )

    if(MKL_LIBRARIES AND MKL_GUIDE)
      set(MKL_LIBRARIES ${MKL_LIBRARIES} mkl_intel mkl_sequential ${MKL_GUIDE} pthread)
    endif()

  endif(${CMAKE_HOST_SYSTEM_PROCESSOR} STREQUAL "x86_64")

endif(CMAKE_MINOR_VERSION GREATER 4)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(MKL DEFAULT_MSG MKL_LIBRARIES)

mark_as_advanced(MKL_LIBRARIES)
mark_as_advanced(MKL_GUIDE)