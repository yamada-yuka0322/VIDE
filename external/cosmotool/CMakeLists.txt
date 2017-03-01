cmake_minimum_required(VERSION 2.6)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}")

project(CosmoToolbox)

include(GetGitRevisionDescription)
include(ExternalProject)
include(FindOpenMP)

get_git_head_revision(HEAD GIT_VER)

option(BUILD_SHARED_LIBS "Build shared libraries." OFF)
option(BUILD_STATIC_LIBS "Build static libraries." ON)
option(ENABLE_OPENMP "Enable OpenMP support." OFF)
option(ENABLE_SHARP "Enable SPHT support." ON)

find_path(NETCDF_INCLUDE_PATH NAMES netcdf.h)
find_path(NETCDFCPP_INCLUDE_PATH NAMES netcdfcpp.h netcdf)
find_path(GSL_INCLUDE_PATH NAMES gsl/gsl_blas.h)

IF(EXISTS ${NETCDFCPP_INCLUDE_PATH}/netcdf)
  SET(FOUND_NETCDF4 1)
  FILE(WRITE ${CMAKE_BINARY_DIR}/src/ctool_netcdf_ver.hpp "#define NETCDFCPP4 1")
ELSE(EXISTS ${NETCDFCPP_INCLUDE_PATH}/netcdf)
  SET(FOUND_NETCDF3 1)
  FILE(WRITE ${CMAKE_BINARY_DIR}/src/ctool_netcdf_ver.hpp "#undef NETCDFCPP4")
ENDIF(EXISTS ${NETCDFCPP_INCLUDE_PATH}/netcdf)

find_library(NETCDF_LIBRARY netcdf)
find_library(NETCDFCPP_LIBRARY NAMES netcdf_c++ netcdf_c++4)
find_library(GSL_LIBRARY gsl)
find_library(GSLCBLAS_LIBRARY gslcblas)
find_library(ZLIB z)
find_library(LDLIB dl)

if (ENABLE_SHARP)
  SET(SHARP_SOURCE ${CMAKE_SOURCE_DIR}/external/sharp)
  SET(DEP_BUILD ${CMAKE_SOURCE_DIR}/external/sharp/auto)
  ExternalProject_Add(sharp
    SOURCE_DIR ${SHARP_SOURCE}
    BUILD_IN_SOURCE 1 
    CONFIGURE_COMMAND ${SHARP_SOURCE}/configure --prefix=${DEP_BUILD}
    BUILD_COMMAND ${CMAKE_MAKE_PROGRAM}
    INSTALL_COMMAND echo "No install"
  )
  SET(SHARP_LIBRARY ${DEP_BUILD}/lib/libsharp.a)
  SET(SHARP_LIBRARIES ${SHARP_LIBRARY})
  SET(SHARP_INCLUDE_PATH ${DEP_BUILD}/include)
endif (ENABLE_SHARP)


find_package(Boost 1.53)
message("INC DIRS: ${Boost_INCLUDE_DIRS}")
mark_as_advanced(Boost_INCLUDE_DIRS Boost_LIBRARIES)


if(HDF5_ROOTDIR)
  SET(ENV{HDF5_ROOT} ${HDF5_ROOTDIR})
endif(HDF5_ROOTDIR)
find_package(HDF5 REQUIRED COMPONENTS HL CXX)

include(FindPkgConfig)

pkg_check_modules(FFTW3 fftw3>=3.3)
pkg_check_modules(FFTW3F fftw3f>=3.3)
pkg_check_modules(EIGEN3 eigen3)

include(FindPackageHandleStandardArgs)
set(NETCDF_FIND_REQUIRED TRUE)
set(GSL_FIND_REQUIRED TRUE)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(NetCDF DEFAULT_MSG NETCDF_LIBRARY NETCDFCPP_LIBRARY NETCDF_INCLUDE_PATH)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(GSL DEFAULT_MSG GSL_LIBRARY GSLCBLAS_LIBRARY GSL_INCLUDE_PATH)


set(GSL_LIBRARIES ${GSL_LIBRARY} ${GSLCBLAS_LIBRARY})

# CPACK Configuration
SET(CPACK_PACKAGE_DESCRIPTION_SUMMARY "A toolbox for impatient cosmologists")
SET(CPACK_PACKAGE_VENDOR "Guilhem Lavaux")
SET(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/LICENCE_CeCILL_V2")
SET(CPACK_PACKAGE_VERSION_MAJOR "1")
SET(CPACK_PACKAGE_VERSION_MINOR "0")
SET(CPACK_PACKAGE_VERSION_PATCH "0-${GIT_VER}")
SET(CPACK_PACKAGE_INSTALL_DIRECTORY "CosmoToolbox-${GalaxExplorer_VERSION_MAJOR}.${GalaxExplorer_VERSION_MINOR}")
SET(CPACK_STRIP_FILES "lib/libCosmoTool.so")
SET(CPACK_SOURCE_IGNORE_FILES 
"/CVS/;/\\\\.git/;/\\\\.svn/;\\\\.swp$;\\\\.#;/#;.*~;cscope.*;/CMakeFiles/;.*\\\\.cmake;Makefile")

add_subdirectory(src)
add_subdirectory(sample)

include(CPack)
