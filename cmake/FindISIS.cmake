# Find ISIS and its dependencies, such as CSPICE, Qt, and patchelf

find_path(ISIS_INCLUDE_DIR
  NAME   Isis.h
  HINTS  "${ISIS_ROOT}/inc"
         "$ENV{ISISROOT}/inc"
         "$ENV{HOME}/projects/isis/inc"
         "$ENV{HOME}/packages/isis/inc"
)

if(ISIS_INCLUDE_DIR)
  message(STATUS "Found the ISIS include directory: ${ISIS_INCLUDE_DIR}")
else()
  message(SEND_ERROR "Could not find the ISIS include directory: ${ISIS_INCLUDE_DIR}")
endif()

if ("${ISIS_ROOT}" STREQUAL "")
  get_filename_component(ISIS_ROOT ${ISIS_INCLUDE_DIR} PATH)
endif()
message(STATUS "ISIS_ROOT is ${ISIS_ROOT}")

# Searching for each library that makes up a component
foreach(COMPONENT ${ISIS_FIND_COMPONENTS})
  string(TOUPPER ${COMPONENT} UPPERCOMPONENT)
  set( ISIS_${UPPERCOMPONENT}_LIBRARY "ISIS_${UPPERCOMPONENT}_LIBRARY-NOTFOUND")

  find_library(ISIS_${UPPERCOMPONENT}_LIBRARY
    NAMES ${COMPONENT}
    HINTS ${ISIS_ROOT}
    PATH_SUFFIXES lib lib64
    )
  
  mark_as_advanced( ISIS_${UPPERCOMPONENT}_LIBRARY )

  if(ISIS_${UPPERCOMPONENT}_LIBRARY)
    set(ISIS_${UPPERCOMPONENT}_FOUND TRUE CACHE INTERNAL "If the ISIS ${UPPERCOMPONENT} library was found")
  endif()

endforeach(COMPONENT)

set(ISIS_INCLUDE_DIRS ${ISIS_INCLUDE_DIR})
set(ISIS_3RD "${ISIS_ROOT}/3rdParty/lib")

if(ISIS_INCLUDE_DIR)
  set(ISIS_FOUND TRUE)
else(ISIS_INCLUDE_DIR)
  set(ISIS_FOUND FALSE)
endif(ISIS_INCLUDE_DIR)

# Summarize
if(ISIS_FOUND)
  message(STATUS "Found the following ISIS libraries:")
  foreach( COMPONENT ${ISIS_FIND_COMPONENTS})
    string(TOUPPER ${COMPONENT} UPPERCOMPONENT )
    if ( ISIS_${UPPERCOMPONENT}_FOUND )
      message(STATUS "  ${ISIS_${UPPERCOMPONENT}_LIBRARY}")
      set(ISIS_LIBRARIES ${ISIS_LIBRARIES} ${ISIS_${UPPERCOMPONENT}_LIBRARY})
    else ( ISIS_${UPPERCOMPONENT}_FOUND )
      message(SEND_ERROR "Unable to find ${COMPONENT}")
    endif()
  endforeach()
else(ISIS_FOUND)
  message(SEND_ERROR "Unable to find requested ISIS libraries")
endif(ISIS_FOUND)

# Search for the cspice dir
find_path(CSPICE_INCLUDE_DIR
  NAMES  SpiceCK.h
  HINTS  "${CSPICE_ROOT}/naif"
         "${ISIS_ROOT}/naif"
         "$ENV{HOME}/packages/cspice/naif"
         "$ENV{HOME}/projects/cspice/naif"
)
if(CSPICE_INCLUDE_DIR)
  message(STATUS "Found the CSPICE include directory: ${CSPICE_INCLUDE_DIR}")
else()
  message(SEND_ERROR "Could not find the CSPICE include directory: ${CSPICE_INCLUDE_DIR}")
endif()

if ("${CSPICE_ROOT}" STREQUAL "")
  get_filename_component(CSPICE_ROOT ${CSPICE_INCLUDE_DIR} PATH)
endif()
message(STATUS "CSPICE_ROOT is ${CSPICE_ROOT}")

if (NOT APPLE)

  # Search for patchelf
  get_filename_component(PATCHELF_GUESS_DIR "${PATCHELF}" PATH)
  find_path(PATCHELF_DIR
    NAMES  patchelf
    HINTS  "${PATCHELF_GUESS_DIR}"
           "${ISIS_ROOT}/bin"
    )
  set(PATCHELF "${PATCHELF_DIR}/patchelf")
  if (EXISTS ${PATCHELF})
    message(STATUS "Found patchelf: ${PATCHELF}")
  else()
    message(SEND_ERROR "Invalid path to patchelf: ${PATCHELF}")
  endif()
  
  # Search for Qt
  find_path(QT_DIR
    NAMES  include/Qt
    HINTS  "${QT_ROOT}"
           "${ISIS_ROOT}"
    )
  set(QT_ROOT "${QT_DIR}")
  if (EXISTS "${QT_ROOT}/")
    message(STATUS "Found Qt in: ${QT_ROOT}")
  else()
    message(SEND_ERROR "Invalid path to Qt: ${QT_ROOT}")
  endif()
  
endif()
  
if (APPLE)

    if (EXISTS "${ISIS_3RD}/libcamd.dylib")
        message(STATUS "Found ${ISIS_3RD}/libcamd.dylib")
        set( CHOLMODLIB_DEPS "-lcamd")
    else()
        set( CHOLMODLIB_DEPS "")
    endif()
        
    set(COMPILE_FLAGS
        "-I${ISIS_INCLUDE_DIR} -I${CSPICE_ROOT} -I${ISIS_3RD}/QtCore.framework/Headers -I${ISIS_3RD}/QtAssistant.framework/Headers -I${ISIS_3RD}/QtGui.framework/Headers -I${ISIS_3RD}/QtNetwork.framework/Headers -I${ISIS_3RD}/QtOpenGL.framework/Headers -I${ISIS_3RD}/QtScript.framework/Headers -I${ISIS_3RD}/QtScriptTools.framework/Headers -I${ISIS_3RD}/QtSql.framework/Headers -I${ISIS_3RD}/QtSvg.framework/Headers -I${ISIS_3RD}/QtTest.framework/Headers -I${ISIS_3RD}/QtWebKit.framework/Headers -I${ISIS_3RD}/QtXml.framework/Headers -I${ISIS_3RD}/QtXmlPatterns.framework/Headers -F${ISIS_3RD}  -I${ISIS_3RD}/qwt.framework/Headers  -I${ISIS_INCLUDE_DIR}/xercesc  -I${ISIS_INCLUDE_DIR}/tiff -I${ISIS_INCLUDE_DIR}/tnt  -I${ISIS_INCLUDE_DIR}/jama -I${ISIS_INCLUDE_DIR}/geos -I${ISIS_INCLUDE_DIR}/gsl -I${ISIS_INCLUDE_DIR}/google -I${ISIS_INCLUDE_DIR}/superlu -Wall -ansi -arch x86_64 -Xarch_x86_64 -mmacosx-version-min=10.6 -DISIS_LITTLE_ENDIAN=1 -D_LARGEFILE64_SOURCE -D_LARGEFILE_SOURCE -fPIC -DGMM_USES_SUPERLU -O1 -DQT_NO_DEBUG -DQT_GUI_LIB -DQT_CORE_LIB -DQT_SHARED -DQT_GUI_LIB -DQT_CORE_LIB -DQT_SHARED -DENABLEJP2K=0")
    set(LINK_FLAGS
        "-L${ISIS_ROOT}/lib -L${ISIS_3RD} -F${ISIS_3RD} -headerpad_max_install_names -arch x86_64 -Xarch_x86_64 -mmacosx-version-min=10.6 -bind_at_load -Wl,-rpath,${ISIS_ROOT} -Wl,-w -Wl,-rpath,@loader_path/.. -lz -lm -framework ApplicationServices -framework QtXmlPatterns -framework QtXml -framework QtNetwork -framework QtSql -framework QtGui -framework QtCore -framework QtSvg -framework QtTest -framework QtWebKit -framework QtOpenGL -framework qwt -lxerces-c -ltiff -lcspice -lgeos -lgsl -lgslcblas -lprotobuf -lkdu_a63R -lcholmod -lamd -lcolamd ${CHOLMODLIB_DEPS} -framework Accelerate -lsuperlu -lblas")
else()
    set(COMPILE_FLAGS "-I${ISIS_ROOT}/inc -I${CSPICE_ROOT} -I${QT_ROOT}/include -I${QT_ROOT}/include/Qt -I${QT_ROOT}/include/QtCore -I${QT_ROOT}/include/QtAssistant -I${QT_ROOT}/include/QtGui -I${QT_ROOT}/include/QtNetwork -I${QT_ROOT}/include/QtOpenGL -I${QT_ROOT}/include/QtScript -I${QT_ROOT}/include/QtScriptTools -I${QT_ROOT}/include/QtSql -I${QT_ROOT}/include/QtSvg -I${QT_ROOT}/include/QtTest -I${QT_ROOT}/include/QtWebKit -I${QT_ROOT}/include/QtXml -I${QT_ROOT}/include/QtXmlPatterns -I${ISIS_ROOT}/3rdParty/include/qwt/qwt6.0.1 -I${ISIS_ROOT}/3rdParty/include/xercesc/xercesc-3.1.1 -I${ISIS_ROOT}/3rdParty/include/tiff/tiff-4.0.1 -I${ISIS_ROOT}/3rdParty/include/naif/cspice64 -I${ISIS_ROOT}/3rdParty/include/tnt/tnt126 -I${ISIS_ROOT}/3rdParty/include/tnt/tnt126/tnt -I${ISIS_ROOT}/3rdParty/include/jama/jama125 -I${ISIS_ROOT}/3rdParty/include/geos/geos3.3.2 -I${ISIS_ROOT}/3rdParty/include/gmm/gmm-4.1 -I${ISIS_ROOT}/3rdParty/include/gmm/gmm-4.1/gmm -I${ISIS_ROOT}/3rdParty/include/google-protobuf/protobuf2.4.1 -I${ISIS_ROOT}/3rdParty/include/boost/boost1.48.0 -I${ISIS_ROOT}/3rdParty/include/kakadu/v6_3-00967N -I${ISIS_ROOT}/3rdParty/include/CHOLMOD/CHOLMOD1.7.3 -I${ISIS_ROOT}/3rdParty/include/superlu/superlu4.3 -Wall -ansi -DISIS_LITTLE_ENDIAN=1 -fPIC -DGMM_USES_SUPERLU -O1 -DENABLEJP2K=0")
  set(LINK_FLAGS
"-L${ISIS_ROOT}/lib -L${ISIS_ROOT}/3rdParty/lib -Wl,-E -Xlinker -z -Xlinker origin -Xlinker -rpath -Xlinker ${ISIS_ROOT}/lib:${ISIS_ROOT}/3rdParty/lib -lisis3.4.6 -pthread -lQtXmlPatterns -lQtXml -lQtNetwork -lQtSql -lQtCore -lQtSvg -lQtTest -lQtOpenGL -lQtDBus -lQtWebKit -lqwt -lxerces-c -ltiff -lcspice -lgeos-3.3.2 -lgsl -lgslcblas -lX11 -lprotobuf -lkdu_a63R -lcholmod -lamd -lcolamd -llapack -lsuperlu_4.3 -lblas -lgfortran")
endif()
    
macro(add_isis_executable name)
  add_executable(${name} ${ARGN})
  
  # Add ISIS libs
  foreach(lib ${ISIS_LIBRARIES})
    target_link_libraries( ${name} ${lib} )
  endforeach(lib)

  # Add other libs
  foreach(lib ${PHOTK_USED_LIBS})
    if( NOT ${lib} STREQUAL "optimized" AND
        NOT ${lib} STREQUAL "debug" )
        target_link_libraries( ${name} ${lib} )
    endif()
  endforeach(lib)
    
  set_target_properties( ${name}
    PROPERTIES
    COMPILE_FLAGS ${COMPILE_FLAGS}
    LINK_FLAGS ${LINK_FLAGS}
  )
  
endmacro(add_isis_executable name)

macro(edit_tool_paths name)

 # Fix the paths to dynamic libraries. This is needed for ISIS.

  set(relocdir ${ISIS_ROOT}/lib:${ISIS_ROOT}/3rdParty/lib/:${ISIS_ROOT}/3rdParty:${ISIS_ROOT}:${PHOTK_BINARY_DIR}/src/photk )
  foreach(lib ${PHOTK_USED_LIBS})
    get_filename_component(lib_dir ${lib} PATH)
    set(relocdir ${relocdir}:${lib_dir})
  endforeach(lib)

  if(APPLE)
    add_custom_command(
        TARGET ${name}
        POST_BUILD
        COMMAND ${ISIS_ROOT}/scripts/SetRunTimePath --bins --libmap=${ISIS_ROOT}/scripts/isis_bins_paths.lis --liblog=${ISIS_ROOT}/scripts/DarwinLibs.lis --relocdir=${relocdir} --errlog=DarwinErrors.lis ${name}
        COMMENT "Editing the run-time paths"
        VERBATIM
    )
  else()
    add_custom_command(
        TARGET ${name}
        POST_BUILD
        COMMAND ${PATCHELF} --set-rpath ${relocdir} ${name}
        COMMENT "Editing the run-time paths"
        VERBATIM
    )
  endif()
endmacro(edit_tool_paths)
  
macro(add_isis_tool name)
  add_isis_executable(${name} ${ARGN})
  if( PHOTK_BUILD_TOOLS )
    message( "Installing: " ${name} )
    install(TARGETS ${name} RUNTIME DESTINATION bin)
  endif()
  set_target_properties(${name} PROPERTIES FOLDER "tools")
  edit_tool_paths(${name})
  
endmacro(add_isis_tool name)

