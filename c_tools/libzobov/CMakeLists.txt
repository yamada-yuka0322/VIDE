
add_library(zobovTool loadZobov.cpp particleInfo.cpp contour_pixels.cpp)

set_target_properties(zobovTool PROPERTIES COMPILE_FLAGS ${OpenMP_CXX_FLAGS} LINK_FLAGS ${OpenMP_CXX_FLAGS})
target_link_libraries(zobovTool ${COSMOTOOL_LIBRARY} ${GSL_LIBRARIES} ${NETCDF_LIBRARIES} ${HDF5_CXX_LIBRARIES} ${HDF5_HL_LIBRARIES} ${HDF5_LIBRARIES} ${DL_LIBRARY})
