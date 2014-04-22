
add_library(zobovTool loadZobov.cpp particleInfo.cpp contour_pixels.cpp)

set_target_properties(zobovTool PROPERTIES COMPILE_FLAGS ${OpenMP_CXX_FLAGS} LINK_FLAGS ${OpenMP_CXX_FLAGS})

SET(ZOB_LIBS zobovTool ${COSMOTOOL_LIBRARY} ${GSL_LIBRARY} ${GSLCBLAS_LIBRARY} ${NETCDFCPP_LIBRARY} ${NETCDF_LIBRARY})
SET(ZOB_LIBS PARENT ${ZOB_LIBS})

