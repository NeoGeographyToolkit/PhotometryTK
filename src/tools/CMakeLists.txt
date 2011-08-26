include_directories(${VISIONWORKBENCH_INCLUDE_DIRS})
include_directories(${Boost_INCLUDE_DIRS})

# Reusable Vision Workbench Patterns
set(VISIONWORKBENCH_BASE_LIBRARIES
  ${VISIONWORKBENCH_CORE_LIBRARY}
  ${VISIONWORKBENCH_MATH_LIBRARY}
  ${VISIONWORKBENCH_IMAGE_LIBRARY}
  ${VISIONWORKBENCH_FILEIO_LIBRARY}
  )
set(VISIONWORKBENCH_BASE_PLATE_LIBRARIES
  ${VISIONWORKBENCH_BASE_LIBRARIES}
  ${VISIONWORKBENCH_PLATE_LIBRARY}
  )

# Executables
set(PHOTK_USED_LIBS ${VISIONWORKBENCH_BASE_PLATE_LIBRARIES} ${Boost_PROGRAM_OPTIONS_LIBRARY} ${Boost_THREAD_LIBRARY})
#add_photk_tool(mipmap
#  mipmap.cc)

set(PHOTK_USED_LIBS PhotometryTK ${VISIONWORKBENCH_BASE_PLATE_LIBRARIES} ${VISIONWORKBENCH_CARTOGRAPHY_LIBRARY} ${PROTOBUF_LIBRARIES} ${Boost_LIBRARIES})
add_photk_tool(phodrg2plate
  phodrg2plate.cc)
add_photk_tool(phoinitet
  phoinitet.cc)
add_photk_tool(phoinitfile
  phoinitfile.cc)
add_photk_tool(phoitalbedo
  phoitalbedo.cc)
add_photk_tool(phoitdem
  phoitdem.cc)
add_photk_tool(phoiterror
  phoiterror.cc)
add_photk_tool(phoitnorm
  phoitnorm.cc)
add_photk_tool(phoitreflect
  phoitreflect.cc)
add_photk_tool(phoittime
  phoittime.cc)

set(PHOTK_USED_LIBS PhotometryTK ${PROTOBUF_LIBRARIES} ${VISIONWORKBENCH_BASE_PLATE_LIBRARIES} ${Boost_LIBRARIES})
add_photk_tool(ptk_server
  ptk_server.cc)

set(PHOTK_USED_LIBS ${VISIONWORKBENCH_BASE_LIBRARIES} ${VISIONWORKBENCH_CARTOGRAPHY_LIBRARY} ${Boost_LIBRARIES} )
add_photk_tool(shadow_mask
  shadow_mask.cc)

# Scripts
install( FILES phoinit.py phosolve.py pmipmap.py
  PERMISSIONS GROUP_EXECUTE OWNER_EXECUTE OWNER_READ GROUP_READ
  DESTINATION bin )