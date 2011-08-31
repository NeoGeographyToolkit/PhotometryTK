include_directories(${VISIONWORKBENCH_INCLUDE_DIRS})
include_directories(${Boost_INCLUDE_DIRS})

include_directories(${CMAKE_CURRENT_BINARY_DIR}/..)
PROTOBUF_GENERATE_CPP(
  PROTO_SRCS
  PROTO_HDRS
  ProjectFile.proto
  )

add_library(PhotometryTK
  Common.cc
  PlateCommon.cc
  ProjectService.cc
  RemoteProjectFile.cc
  ${PROTO_SRCS}
  )

target_link_libraries(PhotometryTK
  ${PROTOBUF_LIBRARIES}
  ${Boost_SYSTEM_LIBRARY}
  ${Boost_FILESYSTEM_LIBRARY}
  ${Boost_THREAD_LIBRARY}
  ${Boost_PROGRAM_OPTIONS_LIBRARY}
  ${VISIONWORKBENCH_CORE_LIBRARY}
  ${VISIONWORKBENCH_MATH_LIBRARY}
  ${VISIONWORKBENCH_IMAGE_LIBRARY}
  ${VISIONWORKBENCH_FILEIO_LIBRARY}
  ${VISIONWORKBENCH_PLATE_LIBRARY}
  )

install(FILES
  ${PROTO_HDRS}
  AlbedoAccumulators.h
  Common.h
  ErrorAccumulators.h
  ImageAccumulators.h
  Macros.h
  PlateCommon.h
  ProjectFileIO.h
  ProjectService.h
  RecursiveBBoxAccumulator.h
  RemoteProjectFile.h
  TimeAccumulators.h
  DESTINATION include/photk
  )

install(TARGETS
  PhotometryTK
  DESTINATION lib
  )