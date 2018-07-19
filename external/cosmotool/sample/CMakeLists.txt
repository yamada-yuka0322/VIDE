SET(tolink  ${CosmoTool_local} ${CosmoTool_LIBS} ${GSL_LIBRARIES} ${DL_LIBRARY})
include_directories(${CMAKE_SOURCE_DIR}/src)
include_directories(${FFTW3_INCLUDE_DIRS} ${EIGEN3_INCLUDE_DIRS} ${GSL_INCLUDE_PATH})
if(YORICK_SUPPORT)
  include_directories(${NETCDF_INCLUDE_PATH})
endif(YORICK_SUPPORT)
include_directories(${CMAKE_SOURCE_DIR}/sample)

IF(SHARP_INCLUDE_PATH)
  include_directories(BEFORE ${SHARP_INCLUDE_PATH})
ENDIF(SHARP_INCLUDE_PATH)

add_executable(testBQueue testBQueue.cpp)
target_link_libraries(testBQueue ${tolink})

add_executable(testInterpolate testInterpolate.cpp)
target_link_libraries(testInterpolate ${tolink})

add_executable(testkd testkd.cpp)
target_link_libraries(testkd ${tolink})

add_executable(testkd2 testkd2.cpp)
target_link_libraries(testkd2 ${tolink})

add_executable(testkd3 testkd3.cpp)
target_link_libraries(testkd3 ${tolink})

add_executable(testDelaunay testDelaunay.cpp)
target_link_libraries(testDelaunay ${tolink})

add_executable(testNewton testNewton.cpp)
target_link_libraries(testNewton ${tolink})

add_executable(testPool testPool.cpp)
target_link_libraries(testPool ${tolink})

if (HDF5_FOUND)
    include_directories(${HDF5_INCLUDE_DIR})

    add_executable(testReadFlash testReadFlash.cpp)
    target_link_libraries(testReadFlash ${tolink})
endif (HDF5_FOUND)


add_executable(testEskow testEskow.cpp)
target_link_libraries(testEskow ${tolink})

add_executable(testAlgo testAlgo.cpp)
target_link_libraries(testAlgo ${tolink})

add_executable(testBSP testBSP.cpp)
target_link_libraries(testBSP ${tolink})

if (FFTW3_FOUND AND FFTW3F_FOUND AND EIGEN3_FOUND)
  IF (FFTW3F_LIBRARY_FULL)
    SET(FFTW3_LIB ${FFTW3F_LIBRARY_FULL})
  ENDIF (FFTW3F_LIBRARY_FULL)
  IF (FFTW3_LIBRARY_FULL)
    SET(FFTW3_LIB ${FFTW3_LIB} ${FFTW3_LIBRARY_FULL})
  ENDIF (FFTW3_LIBRARY_FULL)
  IF(YORICK_SUPPORT)
    add_executable(test_fft_calls test_fft_calls.cpp)
    target_link_libraries(test_fft_calls ${tolink} ${FFTW3_LIB})
  ENDIF(YORICK_SUPPORT)
endif (FFTW3_FOUND AND FFTW3F_FOUND AND  EIGEN3_FOUND)

if (ENABLE_SHARP AND SHARP_LIBRARY AND SHARP_INCLUDE_PATH AND EIGEN3_FOUND)
  include_directories(${SHARP_INCLUDE_PATH})
  add_executable(test_healpix_calls test_healpix_calls.cpp)
  target_link_libraries(test_healpix_calls ${tolink} ${SHARP_LIBRARIES})
  if (ENABLE_OPENMP)
    set_target_properties(test_healpix_calls PROPERTIES COMPILE_FLAGS "${OpenMP_CXX_FLAGS}" LINK_FLAGS "${OpenMP_CXX_FLAGS}")
  endif()
  add_dependencies(test_healpix_calls sharp)
endif (ENABLE_SHARP AND SHARP_LIBRARY AND SHARP_INCLUDE_PATH AND EIGEN3_FOUND)

add_executable(test_cosmopower test_cosmopower.cpp)
target_link_libraries(test_cosmopower ${tolink})

if (Boost_FOUND)
  include_directories(${Boost_INCLUDE_DIRS})

  IF (YORICK_SUPPORT)
    add_executable(testSmooth testSmooth.cpp)
    target_link_libraries(testSmooth ${tolink})

  ENDIF (YORICK_SUPPORT)
  if (HDF5_FOUND)
    add_executable(simple3DFilter simple3DFilter.cpp)
    target_link_libraries(simple3DFilter ${tolink})

    add_executable(simpleDistanceFilter simpleDistanceFilter.cpp)
    target_link_libraries(simpleDistanceFilter ${tolink})

    install(TARGETS simple3DFilter simpleDistanceFilter
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib)

    add_executable(gadgetToArray gadgetToArray.cpp)
    target_link_libraries(gadgetToArray ${tolink})
    
    add_executable(testHDF5 testHDF5.cpp)
    target_link_libraries(testHDF5 ${tolink})

    add_executable(graficToDensity graficToDensity.cpp)
    target_link_libraries(graficToDensity ${tolink})
  endif()
endif (Boost_FOUND)

IF (ENABLE_OPENMP AND YORICK_SUPPORT)
  add_executable(gadgetToDensity gadgetToDensity.cpp)
  target_link_libraries(gadgetToDensity ${tolink})
ENDIF (ENABLE_OPENMP AND YORICK_SUPPORT)

