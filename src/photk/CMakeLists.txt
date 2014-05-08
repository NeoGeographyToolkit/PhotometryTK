include_directories(${VISIONWORKBENCH_INCLUDE_DIRS})
include_directories(${Boost_INCLUDE_DIRS})
include_directories(${PROTOBUF_INCLUDE_DIRS})

include_directories(${CMAKE_CURRENT_BINARY_DIR}/..)

add_library(PhotometryTK
        Albedo.cc    Index.cc    Reconstruct.cc       Shadow.cc            Weights.cc
        Camres.cc    Misc.cc     ReconstructError.cc  Shape.cc
        Exposure.cc  Outlier.cc  Reflectance.cc       ShapeFromShading.cc
  )

target_link_libraries(PhotometryTK
  ${Boost_SYSTEM_LIBRARY}
  ${Boost_FILESYSTEM_LIBRARY}
  ${Boost_THREAD_LIBRARY}
  ${Boost_PROGRAM_OPTIONS_LIBRARY}
  ${VISIONWORKBENCH_CORE_LIBRARY}
  ${VISIONWORKBENCH_MATH_LIBRARY}
  ${VISIONWORKBENCH_IMAGE_LIBRARY}
  ${VISIONWORKBENCH_FILEIO_LIBRARY}
  ${VISIONWORKBENCH_CARTOGRAPHY_LIBRARY}
  ${StereoPipeline_LIBRARIES}
  )

install(FILES
        Albedo.h    Index.h    Reconstruct.h       Shadow.h            Weights.h
        Camres.h    Misc.h     ReconstructError.h  Shape.h
        Exposure.h  Outlier.h  Reflectance.h       ShapeFromShading.h
          DESTINATION include/photk
  )

install(TARGETS
  PhotometryTK
  DESTINATION lib
  )
