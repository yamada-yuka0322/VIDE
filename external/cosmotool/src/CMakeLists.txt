
	SET(CosmoTool_SRCS
		fortran.cpp
		interpolate.cpp
		load_data.cpp
		loadRamses.cpp
		powerSpectrum.cpp
		miniargs.cpp
		growthFactor.cpp
		cosmopower.cpp
	  cic.cpp
	)

	IF (Boost_FOUND)
	  include_directories(${Boost_INCLUDE_DIRS})
  SET(CosmoTool_SRCS
    ${CosmoTool_SRCS}
    loadGadget.cpp
  )
ENDIF (Boost_FOUND)

IF (ENABLE_OPENMP)
ENDIF (ENABLE_OPENMP)

IF (YORICK_SUPPORT)
  IF(FOUND_NETCDF3)
	SET(CosmoTool_SRCS ${CosmoTool_SRCS} yorick_nc3.cpp)
  ELSE(FOUND_NETCDF3)
  IF(FOUND_NETCDF4)
	SET(CosmoTool_SRCS ${CosmoTool_SRCS} yorick_nc4.cpp)
  ENDIF(FOUND_NETCDF4)
  ENDIF(FOUND_NETCDF3)
ENDIF(YORICK_SUPPORT)


if (HDF5_FOUND)
  add_library(CosmoHDF5 OBJECT
    h5_readFlash.cpp
    loadFlash.cpp
    )
  target_include_directories(CosmoHDF5 BEFORE PRIVATE ${HDF5_INCLUDE_DIR})
else(HDF5_FOUND)
  add_library(CosmoHdf5  OBJECT
    loadFlash_dummy.cpp
  )
endif (HDF5_FOUND)

SET(CosmoTool_SRCS ${CosmoTool_SRCS}
	bqueue.hpp
	config.hpp
	dinterpolate.hpp
	field.hpp
	fixArray.hpp
	fortran.hpp
	interpolate3d.hpp
	interpolate.hpp
	kdtree_leaf.hpp
	load_data.hpp
	loadGadget.hpp
	loadRamses.hpp
	loadSimu.hpp
	miniargs.hpp
	mykdtree.hpp
	octTree.hpp
	powerSpectrum.hpp
	sparseGrid.hpp
	sphSmooth.hpp
	yorick.hpp
	growthFactor.hpp
)

include_directories(${GSL_INCLUDE_PATH} ${CMAKE_BINARY_DIR}/src)

set(CosmoTool_LIBS ${GSL_LIBRARIES})
if(YORICK_SUPPORT)
  set(CosmoTool_LIBS ${CosmoTool_LIBS} ${NETCDFCPP_LIBRARY} ${NETCDF_LIBRARY} ${GSL_LIBRARIES})
  include_directories(${NETCDF_INCLUDE_PATH} ${NETCDFCPP_INCLUDE_PATH})
endif(YORICK_SUPPORT)

if (HDF5_FOUND)
  set(CosmoTool_LIBS ${CosmoTool_LIBS} ${HDF5_CXX_LIBRARIES} ${HDF5_HL_LIBRARIES} ${HDF5_LIBRARIES} ${ZLIB})
  message(STATUS "CosmoTool_LIBS: ${CosmoTool_LIBS}")
endif (HDF5_FOUND)

set(CosmoTool_LIBS ${CosmoTool_LIBS} PARENT_SCOPE)

if (BUILD_SHARED_LIBS)
  add_library(CosmoTool SHARED ${CosmoTool_SRCS} $<TARGET_OBJECTS:CosmoHDF5>)
  target_link_libraries(CosmoTool ${CosmoTool_LIBS})
  if (BUILD_STATIC_LIBS)
     add_library(CosmoTool_static STATIC ${CosmoTool_SRCS} $<TARGET_OBJECTS:CosmoHDF5>)
     target_link_libraries(CosmoTool_static)
     set_property(TARGET CosmoTool_static PROPERTY POSITION_INDEPENDENT_CODE ON)
  endif(BUILD_STATIC_LIBS)
else (BUILD_SHARED_LIBS)
  add_library(CosmoTool STATIC ${CosmoTool_SRCS} $<TARGET_OBJECTS:CosmoHDF5>)
  target_link_libraries(CosmoTool)
  set_property(TARGET CosmoTool PROPERTY POSITION_INDEPENDENT_CODE ON)
endif (BUILD_SHARED_LIBS)

install(TARGETS CosmoTool 
	LIBRARY DESTINATION lib
	ARCHIVE DESTINATION lib)

if (BUILD_SHARED_LIBS)
  install(TARGETS CosmoTool_static 
  	  LIBRARY DESTINATION lib
	  ARCHIVE DESTINATION lib)
endif (BUILD_SHARED_LIBS)

install(DIRECTORY . DESTINATION include/CosmoTool
	FILES_MATCHING PATTERN "*.hpp")
install(DIRECTORY ${OMPTL_BUILD_DIR}
     DESTINATION include/CosmoTool)
install(DIRECTORY . DESTINATION include/CosmoTool
	FILES_MATCHING PATTERN "*.tcc")
