include_directories(${VISIONWORKBENCH_INCLUDE_DIRS})
include_directories(${Boost_INCLUDE_DIRS})
 
set(VISIONWORKBENCH_BASE_LIBRARIES
  ${VISIONWORKBENCH_CORE_LIBRARY}
  ${VISIONWORKBENCH_MATH_LIBRARY}
  ${VISIONWORKBENCH_IMAGE_LIBRARY}
  ${VISIONWORKBENCH_FILEIO_LIBRARY}
  )

set(PHOTK_USED_LIBS PhotometryTK ${VISIONWORKBENCH_CARTOGRAPHY_LIBRARY} ${Boost_LIBRARIES})

# photk_add_test(TestAlbedoAccumulators.cxx)
# photk_add_test(TestBBoxIterator.cxx)
# photk_add_test(TestCommon.cxx)
# photk_add_test(TestErrorAccumulators.cxx)
photk_add_test(TestFullIteration.cxx)
# photk_add_test(TestPlateCommon.cxx)
# photk_add_test(TestProjectFileIO.cxx)
# photk_add_test(TestProjectService.cxx)
# photk_add_test(TestRecursiveBBoxAccumulator.cxx)
# photk_add_test(TestTimeAccumulators.cxx)


