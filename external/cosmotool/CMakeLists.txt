cmake_minimum_required(VERSION 3.6)
set(CMAKE_CXX_STANDARD 11)
list(INSERT CMAKE_MODULE_PATH 0 "${CMAKE_CURRENT_SOURCE_DIR}")

project(CosmoToolbox)

include(GetGitRevisionDescription)
include(ExternalProject)
include(FindPkgConfig)
include(FindPackageHandleStandardArgs)
include(color_msg)

option(BUILD_SHARED_LIBS "Build shared libraries." OFF)
option(BUILD_STATIC_LIBS "Build static libraries." ON)
option(ENABLE_SHARP "Enable SHARP support." ON)
option(INTERNAL_HDF5 "Build internal version of HDF5" OFF)
option(INTERNAL_EIGEN "Build internal version of Eigen" ON)
option(INTERNAL_NETCDF "Build internal version of NETCDF" OFF)
option(INTERNAL_BOOST "Build internal version of BOOST" OFF)
option(USE_GIT_VERSION "Use GIT version to pad version number" OFF)
option(YORICK_SUPPORT "Add for pleasant yorick support" ON)
option(BUILD_PYTHON "Build the python modules" ON)
IF(USE_GIT_VERSION)
  get_git_head_revision(HEAD GIT_VER)
  SET(EXTRA_VERSION "-${GIT_VER}")
ENDIF(USE_GIT_VERSION)

IF(BUILD_PYTHON)
  include(FindPythonInterp)
  set(NumPy_FIND_REQUIRED TRUE)
  include(${CMAKE_SOURCE_DIR}/FindPyLibs.cmake)
  include(${CMAKE_SOURCE_DIR}/FindNumPy.cmake)
ENDIF(BUILD_PYTHON)


#IF(NOT BUILD_SHARED_LIBS AND BUILD_STATIC_LIBS)
#  SET(CosmoTool_local CosmoTool_static)
#ELSE(NOT BUILD_SHARED_LIBS AND BUILD_STATIC_LIBS)
  SET(CosmoTool_local CosmoTool)
#ENDIF(NOT BUILD_SHARED_LIBS AND BUILD_STATIC_LIBS)

MESSAGE(STATUS "Using the target ${CosmoTool_local} to build python module")

include(${CMAKE_SOURCE_DIR}/external/external_build.cmake)

IF(YORICK_SUPPORT)
  IF((EXISTS ${NETCDFCPP_INCLUDE_PATH}/netcdf AND ${NETCDFCPP_LIBRARY} MATCHES "netcdf_c\\+\\+4") OR (INTERNAL_NETCDF))
    SET(FOUND_NETCDF4 1)
    FILE(WRITE ${CMAKE_BINARY_DIR}/src/ctool_netcdf_ver.hpp "#define NETCDFCPP4 1")
  ELSE((EXISTS ${NETCDFCPP_INCLUDE_PATH}/netcdf AND ${NETCDFCPP_LIBRARY} MATCHES "netcdf_c\\+\\+4") OR (INTERNAL_NETCDF))
    SET(FOUND_NETCDF3 1)
    FILE(WRITE ${CMAKE_BINARY_DIR}/src/ctool_netcdf_ver.hpp "#undef NETCDFCPP4")
  ENDIF((EXISTS ${NETCDFCPP_INCLUDE_PATH}/netcdf AND ${NETCDFCPP_LIBRARY} MATCHES "netcdf_c\\+\\+4") OR (INTERNAL_NETCDF))
ENDIF(YORICK_SUPPORT)

find_program(CYTHON cython3 cython)
find_library(ZLIB z)
find_library(DL_LIBRARY dl)
find_library(MATH_LIBRARY m)

set(NETCDF_FIND_REQUIRED ${YORICK_SUPPORT})
set(GSL_FIND_REQUIRED TRUE)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(NetCDF DEFAULT_MSG NETCDF_LIBRARY NETCDFCPP_LIBRARY NETCDF_INCLUDE_PATH)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(GSL DEFAULT_MSG GSL_LIBRARY GSLCBLAS_LIBRARY GSL_INCLUDE_PATH)


# CPACK Configuration
SET(CPACK_PACKAGE_DESCRIPTION_SUMMARY "A toolbox for impatient cosmologists")
SET(CPACK_PACKAGE_VENDOR "Guilhem Lavaux")
SET(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/LICENCE_CeCILL_V2")
SET(CPACK_PACKAGE_VERSION_MAJOR "1")
SET(CPACK_PACKAGE_VERSION_MINOR "0")
SET(CPACK_PACKAGE_VERSION_PATCH "0${EXTRA_VERSION}")
SET(CPACK_PACKAGE_INSTALL_DIRECTORY "CosmoToolbox-${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}")
SET(CPACK_STRIP_FILES "lib/libCosmoTool.so")
SET(CPACK_SOURCE_IGNORE_FILES 
"/CVS/;/\\\\.git/;/\\\\.svn/;\\\\.swp$;\\\\.#;/#;.*~;cscope.*;/CMakeFiles/;.*\\\\.cmake;Makefile")

add_subdirectory(src)
add_subdirectory(sample)
if(BUILD_PYTHON)
  add_subdirectory(python)
endif(BUILD_PYTHON)

include(CPack)
