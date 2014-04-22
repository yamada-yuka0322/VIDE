include_directories(${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/loaders)

IF(SDF_SUPPORT)
  add_definitions(-DSDF_SUPPORT)
ENDIF(SDF_SUPPORT)

SET(generateMock_SRCS generateMock.cpp )
add_genopt(generateMock_SRCS generateMock.ggo generateMock_conf STRUCTNAME generateMock_info)
add_executable(generateMock ${generateMock_SRCS})
target_link_libraries(generateMock simu_loaders ${ZOB_LIBS} ${HDF5HL_CPP_LIBRARY} ${HDF5_CPP_LIBRARY} ${HDF5_LIBRARY} ${LIBSDF_LIBRARY})

add_executable(generateTestMock generateTestMock.cpp)
target_link_libraries(generateTestMock ${ZOB_LIBS})


SET(generateFromCatalog_SRCS generateFromCatalog.cpp)
add_genopt(generateFromCatalog_SRCS generateFromCatalog.ggo generateFromCatalog_conf STRUCTNAME generateFromCatalog_info)
add_executable(generateFromCatalog ${generateFromCatalog_SRCS})
target_link_libraries(generateFromCatalog ${ZOB_LIBS} ${HDF5HL_CPP_LIBRARY} ${HDF5_CPP_LIBRARY} ${HEALPIX_LIBRARIES})

set_target_properties(generateFromCatalog PROPERTIES COMPILE_FLAGS ${OpenMP_CXX_FLAGS} LINK_FLAGS ${OpenMP_CXX_FLAGS})

subdirs(loaders)
