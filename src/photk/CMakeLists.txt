include_directories(${VISIONWORKBENCH_INCLUDE_DIRS})
include_directories(${Boost_INCLUDE_DIRS})

include_directories(${CMAKE_CURRENT_BINARY_DIR}/..)
PROTOBUF_GENERATE_CPP(
  PROTO_SRCS
  PROTO_HDRS
  ProjectFile.proto
  )

add_library(PhotometryTK
  ProjectService.cc
  RemoteProjectFile.cc
  ${PROTO_SRCS}
  )

install(FILES
  ${PROTO_HDRS}
  AlbedoAccumulators.h
  Common.h
  ErrorAccumulators.h
  ImageAccumulators.h
  Macros.h
  ProjectFileIO.h
  ProjectService.h
  RecursiveBBoxAccumulator.h
  RemoteProjectFile.h
  TimeAccumulators.h
  DESTINATION include/photk
  )