SET(CosmoTool_SRCS
	fortran.cpp
	interpolate.cpp
	load_data.cpp
	loadGadget.cpp
	loadRamses.cpp
	octTree.cpp
	powerSpectrum.cpp
	miniargs.cpp
	growthFactor.cpp
	cosmopower.cpp
)

IF(FOUND_NETCDF3)
	SET(CosmoTool_SRCS ${CosmoTool_SRCS} yorick_nc3.cpp)
ELSE(FOUND_NETCDF3)
IF(FOUND_NETCDF4)
	SET(CosmoTool_SRCS ${CosmoTool_SRCS} yorick_nc4.cpp)
ENDIF(FOUND_NETCDF4)
ENDIF(FOUND_NETCDF3)


if (HDF5_FOUND)
  set(CosmoTool_SRCS ${CosmoTool_SRCS}
    h5_readFlash.cpp
    loadFlash.cpp
    )
else(HDF5_FOUND)
  set(CosmoTool_SRCS ${CosmoTool_SRCS}
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

include_directories(${GSL_INCLUDE_PATH} ${NETCDF_INCLUDE_PATH} ${NETCDFCPP_INCLUDE_PATH} ${CMAKE_BINARY_DIR}/src)

set(CosmoTool_LIBS  ${NETCDFCPP_LIBRARY} ${NETCDF_LIBRARY} ${GSL_LIBRARIES})
if (HDF5_FOUND)
  set(CosmoTool_LIBS ${CosmoTool_LIBS} ${HDF5_CXX_LIBRARIES} ${HDF5_LIBRARIES})
  include_directories(${HDF5_INCLUDE_DIRS})
endif (HDF5_FOUND)

set(CosmoTool_LIBS ${CosmoTool_LIBS} PARENT_SCOPE)

if (BUILD_SHARED_LIBS)
  add_library(CosmoTool SHARED ${CosmoTool_SRCS})
  target_link_libraries(CosmoTool ${CosmoTool_LIBS})

  if (BUILD_STATIC_LIBS)
     add_library(CosmoTool_static STATIC ${CosmoTool_SRCS})
  endif(BUILD_STATIC_LIBS)
else (BUILD_SHARED_LIBS)
  add_library(CosmoTool STATIC ${CosmoTool_SRCS})
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
install(DIRECTORY . DESTINATION include/CosmoTool
	FILES_MATCHING PATTERN "*.tcc")
