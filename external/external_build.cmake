include(FindOpenMP)

OPTION(ENABLE_OPENMP "Set to Yes if Healpix and/or you need openMP" OFF)

IF(ENABLE_OPENMP)

  IF (NOT OPENMP_FOUND)
    MESSAGE(ERROR "No known compiler option for enabling OpenMP")
  ENDIF(NOT OPENMP_FOUND)

  SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKED_FLAGS} ${OpenMP_C_FLAGS}")

ENDIF(ENABLE_OPENMP)

SET(BUILD_PREFIX ${CMAKE_BINARY_DIR}/ep_build)

SET(INTERNAL_FFTW OFF)
SET(INTERNAL_GSL ON)
SET(INTERNAL_BOOST ON)
SET(INTERNAL_NETCDF ON)
SET(INTERNAL_HDF5 ON)
SET(INTERNAL_GENGETOPT ON)
SET(INTERNAL_QHULL ON)

IF(INTERNAL_GENGETOPT)
  SET(GENGETOPT_URL "ftp://ftp.gnu.org/gnu/gengetopt/gengetopt-2.22.5.tar.gz" CACHE STRING "URL to download gengetopt from")
ENDIF(INTERNAL_GENGETOPT)

IF(INTERNAL_HDF5)
  SET(HDF5_URL "http://www.hdfgroup.org/ftp/HDF5/current/src/hdf5-1.8.9.tar.gz" CACHE STRING "URL to download HDF5 from")
ENDIF(INTERNAL_HDF5)

IF(INTERNAL_NETCDF)
  SET(NETCDF_URL "http://www.unidata.ucar.edu/downloads/netcdf/ftp/netcdf-4.1.3.tar.gz" CACHE STRING "URL to download NetCDF from")
ENDIF(INTERNAL_NETCDF)

IF(INTERNAL_BOOST)
  SET(BOOST_URL "http://sourceforge.net/projects/boost/files/boost/1.49.0/boost_1_49_0.tar.gz/download" CACHE STRING "URL to download Boost from")
ELSE(INTERNAL_BOOST)
  find_package(Boost 1.49.0 COMPONENTS format spirit phoenix python FATAL_ERROR)
ENDIF(INTERNAL_BOOST)

IF(INTERNAL_GSL)
  SET(GSL_URL "ftp://ftp.gnu.org/gnu/gsl/gsl-1.15.tar.gz" CACHE STRING "URL to download GSL from ")
ENDIF(INTERNAL_GSL)

IF(INTERNAL_QHULL)
  SET(QHULL_URL "http://www.qhull.org/download/qhull-2012.1-src.tgz" CACHE STRING "URL to download QHull from")
ENDIF(INTERNAL_QHULL)


find_library(ZLIB_LIBRARY z)

SET(CONFIGURE_CPP_FLAGS "${EXTRA_CPP_FLAGS}")
SET(CONFIGURE_LD_FLAGS "${EXTRA_LD_FLAGS}")


##################
# Build gengetopt
##################

if (INTERNAL_GENGETOPT)
  SET(GENGETOPT_SOURCE_DIR ${BUILD_PREFIX}/gengetopt-prefix/src/gengetopt)
  SET(GENGETOPT_BIN_DIR ${CMAKE_BINARY_DIR}/ext_build/gengetopt)
  ExternalProject_Add(gengetopt
    PREFIX ${BUILD_PREFIX}/gengetopt-prefix
    URL ${GENGETOPT_URL}
    CONFIGURE_COMMAND ${GENGETOPT_SOURCE_DIR}/configure 
                --prefix=${GENGETOPT_BIN_DIR} 
		CPPFLAGS=${CONFIGURE_CPP_FLAGS}
		LDFLAGS=${CONFIGURE_LD_FLAGS}
		CC=${CMAKE_C_COMPILER} 
		CXX=${CMAKE_CXX_COMPILER}
    BUILD_IN_SOURCE 1
    INSTALL_COMMAND make install
  )
  SET(GENGETOPT ${GENGETOPT_BIN_DIR}/bin/gengetopt)
else(INTERNAL_GENGETOPT)
  find_program(GENGETOPT gengetopt)
endif(INTERNAL_GENGETOPT)

###############
# Build HDF5
###############

if (INTERNAL_HDF5)
  SET(HDF5_SOURCE_DIR ${BUILD_PREFIX}/hdf5-prefix/src/hdf5)
  SET(HDF5_BIN_DIR ${CMAKE_BINARY_DIR}/ext_build/hdf5)
  ExternalProject_Add(hdf5
    PREFIX ${BUILD_PREFIX}/hdf5-prefix
    URL ${HDF5_URL}
    CONFIGURE_COMMAND ${HDF5_SOURCE_DIR}/configure --disable-shared --enable-cxx --with-pic --prefix=${HDF5_BIN_DIR} CPPFLAGS=${CONFIGURE_CPP_FLAGS} CC=${CMAKE_C_COMPILER} CXX=${CMAKE_CXX_COMPILER}
    BUILD_IN_SOURCE 1
    INSTALL_COMMAND make install
  )
  SET(cosmotool_DEPS ${cosmotool_DEPS} hdf5)
  SET(hdf5_built hdf5)
  set(HDF5_LIBRARY ${HDF5_BIN_DIR}/lib/libhdf5.a CACHE STRING "HDF5 lib" FORCE)
  set(HDF5_CPP_LIBRARY ${HDF5_BIN_DIR}/lib/libhdf5_cpp.a CACHE STRING "HDF5 C++ lib" FORCE)
  set(HDF5HL_LIBRARY ${HDF5_BIN_DIR}/lib/libhdf5_hl.a CACHE STRING "HDF5-HL lib" FORCE)
  set(HDF5HL_CPP_LIBRARY ${HDF5_BIN_DIR}/lib/libhdf5_hl_cpp.a CACHE STRING "HDF5-HL C++ lib" FORCE)
  SET(HDF5_INCLUDE_PATH ${HDF5_BIN_DIR}/include CACHE STRING "HDF5 include path" FORCE)
  SET(ENV{HDF5_ROOT} ${HDF5_BIN_DIR})
  SET(HDF5_ROOTDIR ${HDF5_BIN_DIR})
  SET(CONFIGURE_LDFLAGS "${CONFIGURE_LDFLAGS} -L${HDF5_BIN_DIR}/lib")
else(INTERNAL_HDF5)
  find_path(HDF5_INCLUDE_PATH hdf5.h)
  find_library(HDF5_LIBRARY hdf5)
  find_library(HDF5_CPP_LIBRARY hdf5_cpp)
  find_library(HDF5HL_CPP_LIBRARY hdf5_hl_cpp)
  find_library(HDF5HL_LIBRARY hdf5_hl)
endif (INTERNAL_HDF5)
SET(CONFIGURE_CPP_FLAGS "${CONFIGURE_CPP_FLAGS} -I${HDF5_INCLUDE_PATH}")

###############
# Build NetCDF
###############


if (INTERNAL_NETCDF)
  SET(NETCDF_SOURCE_DIR ${BUILD_PREFIX}/netcdf-prefix/src/netcdf)
  SET(NETCDF_BIN_DIR ${CMAKE_BINARY_DIR}/ext_build/netcdf)
  SET(CONFIGURE_CPP_FLAGS "${CONFIGURE_CPP_FLAGS} -I${NETCDF_BIN_DIR}/include")
  SET(CONFIGURE_LDFLAGS "${CONFIGURE_LDFLAGS} -L${NETCDF_BIN_DIR}/lib")
  SET(EXTRA_NC_FLAGS CPPFLAGS=${CONFIGURE_CPP_FLAGS} LDFLAGS=${CONFIGURE_LDFLAGS})
  ExternalProject_Add(netcdf
    DEPENDS ${hdf5_built}
    PREFIX ${BUILD_PREFIX}/netcdf-prefix
    URL ${NETCDF_URL}
    CONFIGURE_COMMAND ${NETCDF_SOURCE_DIR}/configure --prefix=${NETCDF_BIN_DIR} --enable-netcdf-4  --with-pic --disable-shared --disable-dap --disable-cdmremote --disable-rpc --disable-examples ${EXTRA_NC_FLAGS} CC=${CMAKE_C_COMPILER} CXX=${CMAKE_CXX_COMPILER}
    BUILD_IN_SOURCE 1
    INSTALL_COMMAND make install
  )
  SET(CONFIGURE_CPP_LDFLAGS "${CONFIGURE_LDFLAGS}")
  SET(EXTRA_NC_FLAGS CPPFLAGS=${CONFIGURE_CPP_FLAGS} LDFLAGS=${CONFIGURE_CPP_LDFLAGS})
  SET(cosmotool_DEPS ${cosmotool_DEPS} netcdf)
  SET(NETCDF_LIBRARY ${NETCDF_BIN_DIR}/lib/libnetcdf.a CACHE STRING "NetCDF lib" FORCE)
  SET(NETCDFCPP_LIBRARY ${NETCDF_BIN_DIR}/lib/libnetcdf_c++.a CACHE STRING "NetCDF-C++ lib" FORCE)
  SET(NETCDF_INCLUDE_PATH ${NETCDF_BIN_DIR}/include CACHE STRING "NetCDF include" FORCE)
  SET(NETCDFCPP_INCLUDE_PATH ${NETCDF_INCLUDE_PATH} CACHE STRING "NetCDF C++ include path" FORCE)

ELSE(INTERNAL_NETCDF)
  find_library(NETCDF_LIBRARY netcdf)
  find_library(NETCDFCPP_LIBRARY netcdf_c++)
  find_path(NETCDF_INCLUDE_PATH NAMES netcdf.h)
  find_path(NETCDFCPP_INCLUDE_PATH NAMES netcdf)
  SET(CONFIGURE_CPP_FLAGS ${CONFIGURE_CPP_FLAGS} -I${NETCDF_INCLUDE_PATH} -I${NETCDFCPP_INCLUDE_PATH})
endif (INTERNAL_NETCDF)

##################
# Build BOOST
##################

if (INTERNAL_BOOST)
  SET(BOOST_SOURCE_DIR ${BUILD_PREFIX}/boost-prefix/src/boost)
  ExternalProject_Add(boost
    URL ${BOOST_URL}
    PREFIX ${BUILD_PREFIX}/boost-prefix
    CONFIGURE_COMMAND ${BOOST_SOURCE_DIR}/bootstrap.sh --prefix=${CMAKE_BINARY_DIR}/ext_build/boost
    BUILD_IN_SOURCE 1
    BUILD_COMMAND ${BOOST_SOURCE_DIR}/b2 --with-exception --with-python
    INSTALL_COMMAND echo "No install"
  )
  set(Boost_INCLUDE_DIRS ${BOOST_SOURCE_DIR} CACHE STRING "Boost path" FORCE)
  set(Boost_LIBRARIES ${BOOST_SOURCE_DIR}/stage/lib/libboost_python.a)
endif (INTERNAL_BOOST)

##################
# Build GSl
##################

IF(INTERNAL_GSL)
  SET(GSL_SOURCE_DIR ${BUILD_PREFIX}/gsl-prefix/src/gsl)
  ExternalProject_Add(gsl
    URL ${GSL_URL}
    PREFIX ${BUILD_PREFIX}/gsl-prefix
    CONFIGURE_COMMAND ${GSL_SOURCE_DIR}/configure --prefix=${CMAKE_BINARY_DIR}/ext_build/gsl --disable-shared CPPFLAGS=${CONFIGURE_CPP_FLAGS} CC=${CMAKE_C_COMPILER} CXX=${CMAKE_CXX_COMPILER}
    BUILD_IN_SOURCE 1
    BUILD_COMMAND make
    INSTALL_COMMAND make install
  )
  SET(GSL_INTERNAL_LIBS ${CMAKE_BINARY_DIR}/ext_build/gsl/lib)
  SET(GSL_LIBRARY ${GSL_INTERNAL_LIBS}/libgsl.a CACHE STRING "GSL internal path" FORCE)
  SET(GSLCBLAS_LIBRARY ${GSL_INTERNAL_LIBS}/libgslcblas.a CACHE STRING "GSL internal path" FORCE)
  set(GSL_INCLUDE_PATH ${CMAKE_BINARY_DIR}/ext_build/gsl/include CACHE STRING "GSL internal path" FORCE)
  SET(cosmotool_DEPS ${cosmotool_DEPS} gsl)
ELSE(INTERNAL_GSL)
  find_library(GSL_LIBRARY gsl)
  find_library(GSLCBLAS_LIBRARY gslcblas)
  find_path(GSL_INCLUDE_PATH NAMES gsl/gsl_blas.h)
ENDIF(INTERNAL_GSL)

##################
# Build CosmoTool
##################


ExternalProject_Add(cosmotool
  DEPENDS ${cosmotool_DEPS}
  PREFIX ${BUILD_PREFIX}/cosmotool-prefix
  SOURCE_DIR ${CMAKE_SOURCE_DIR}/external/cosmotool
  CMAKE_ARGS 
        -DHDF5_DIR=${HDF5_ROOTDIR}
	-DHDF5_ROOTDIR=${HDF5_ROOTDIR}
	-DCMAKE_INSTALL_PREFIX=${CMAKE_BINARY_DIR}/ext_build/cosmotool
	-DNETCDF_INCLUDE_PATH=${NETCDF_INCLUDE_PATH}
	-DNETCDFCPP_INCLUDE_PATH=${NETCDFCPP_INCLUDE_PATH}
	-DGSL_INCLUDE_PATH=${GSL_INCLUDE_PATH} 
	-DGSL_LIBRARY=${GSL_LIBRARY} 
	-DGSLCBLAS_LIBRARY=${GSLCBLAS_LIBRARY}
	-DNETCDF_LIBRARY=${NETCDF_LIBRARY}
	-DNETCDFCPP_LIBRARY=${NETCDFCPP_LIBRARY}
)
SET(COSMOTOOL_LIBRARY ${CMAKE_BINARY_DIR}/ext_build/cosmotool/lib/libCosmoTool.a)
set(COSMOTOOL_INCLUDE_PATH ${CMAKE_BINARY_DIR}/ext_build/cosmotool/include)

#################
# Build cfitsio
#################
ExternalProject_Add(cfitsio
  SOURCE_DIR ${CMAKE_SOURCE_DIR}/external/cfitsio
  PREFIX ${BUILD_PREFIX}/cfitsio-prefix
  CONFIGURE_COMMAND 
        ${CMAKE_SOURCE_DIR}/external/cfitsio/configure 
	    --prefix=${CMAKE_BINARY_DIR}/ext_build/cfitsio 
	    CPPFLAGS=${CONFIGURE_CPP_FLAGS} 
	    CC=${CMAKE_C_COMPILER} 
	    CXX=${CMAKE_CXX_COMPILER}
  BUILD_COMMAND make
  BUILD_IN_SOURCE 1
  INSTALL_COMMAND make install
)
SET(CFITSIO_LIBRARY ${CMAKE_BINARY_DIR}/ext_build/cfitsio/lib/libcfitsio.a)

#################
# Build Healpix 
#################

ExternalProject_Add(healpix
  DEPENDS cfitsio
  PREFIX ${BUILD_PREFIX}/healpix-prefix
  SOURCE_DIR ${CMAKE_SOURCE_DIR}/external/healpix
  CONFIGURE_COMMAND echo No configure
  BUILD_COMMAND make
      HEALPIX_TARGET=sampler 
      HEALPIX_CC=${CMAKE_C_COMPILER} 
      HEALPIX_CXX=${CMAKE_CXX_COMPILER} 
      HEALPIX_BASE_PATH=${CMAKE_BINARY_DIR} 
      OMP_SUPPORT=${ENABLE_OPENMP} 
      EXTRA_CPPFLAGS=${CONFIGURE_CPP_FLAGS} 
      OMP_FLAGS=${OpenMP_C_FLAGS}
  BUILD_IN_SOURCE 1
  INSTALL_COMMAND ${CMAKE_COMMAND} -DHEALPIX_DIR:STRING=${CMAKE_SOURCE_DIR}/external/healpix -DDEST_DIR:STRING=${CMAKE_BINARY_DIR}/ext_build/healpix -P ${CMAKE_SOURCE_DIR}/external/install_healpix.cmake
)
set(HPIX_LIBPATH ${CMAKE_BINARY_DIR}/ext_build/healpix/lib)
set(HEALPIX_LIBRARY ${HPIX_LIBPATH}/libhealpix_cxx.a)
set(FFTPACK_LIBRARY ${HPIX_LIBPATH}/libfftpack.a)
set(CXXSUPPORT_LIBRARY ${HPIX_LIBPATH}/libcxxsupport.a)
set(PSHT_LIBRARY ${HPIX_LIBPATH}/libpsht.a)
set(CUTILS_LIBRARY ${HPIX_LIBPATH}/libc_utils.a)

SET(HEALPIX_INCLUDE_PATH ${CMAKE_BINARY_DIR}/ext_build/healpix/include)
SET(HEALPIX_LIBRARIES ${HEALPIX_LIBRARY} ${FFTPACK_LIBRARY} ${CXXSUPPORT_LIBRARY}  ${PSHT_LIBRARY} ${CUTILS_LIBRARY} ${CFITSIO_LIBRARY} )
set(GSL_LIBRARIES ${GSL_LIBRARY} ${GSLCBLAS_LIBRARY})
SET(NETCDF_LIBRARIES ${NETCDFCPP_LIBRARY} ${NETCDF_LIBRARY} ${HDF5HL_LIBRARY} ${HDF5_LIBRARY} ${ZLIB_LIBRARY})

###############
# Build QHull
###############
if (INTERNAL_QHULL)
  ExternalProject_Add(qhull
    URL ${QHULL_URL}
    PREFIX ${BUILD_PREFIX}/qhull-prefix
    CMAKE_ARGS 
      -DCMAKE_INSTALL_PREFIX=${CMAKE_BINARY_DIR}/ext_build/qhull
      -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
      -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
      
  )
  SET(QHULL_DIR ${CMAKE_BINARY_DIR}/ext_build/qhull)
  SET(QHULL_LIBRARY ${QHULL_DIR}/lib/libqhullstatic_p.a)
  SET(QHULL_CPP_LIBRARY ${QHULL_DIR}/lib/libqhullcpp.a)
  SET(QHULL_INCLUDE_PATH ${QHULL_DIR}/include)

  add_definitions(-Dqh_QHpointer)

else(INTERNAL_QHULL)
endif(INTERNAL_QHULL)

SET(QHULL_LIBRARIES ${QHULL_CPP_LIBRARY} ${QHULL_LIBRARY} )

include_directories(${CMAKE_BINARY_DIR}/src 
                    ${NETCDF_INCLUDE_PATH} ${GSL_INCLUDE_PATH} 
                    ${HDF5_INCLUDE_PATH} ${COSMOTOOL_INCLUDE_PATH} 
                    ${Boost_INCLUDE_DIRS}
		    ${QHULL_INCLUDE_PATH})

