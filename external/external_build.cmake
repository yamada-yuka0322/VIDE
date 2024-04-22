include(FindOpenMP)

OPTION(ENABLE_OPENMP "Set to Yes if Healpix and/or you need openMP" OFF)
OPTION(SDF_SUPPORT "Set to Yes to activate support for SDF" ON)

IF(ENABLE_OPENMP)

  IF (NOT OPENMP_FOUND)
    MESSAGE(WARNING "No known compiler option for enabling OpenMP")
    SET(ENABLE_OPENMP FALSE)
  ENDIF(NOT OPENMP_FOUND)

ENDIF(ENABLE_OPENMP)

SET(BUILD_PREFIX ${CMAKE_BINARY_DIR}/ep_build)
SET(EXT_INSTALL ${CMAKE_BINARY_DIR}/ext_install)

SET(INTERNAL_FFTW OFF)
OPTION(INTERNAL_GSL "Use internal GSL" ON)
OPTION(INTERNAL_BOOST "Use internal Boost" ON)
OPTION(INTERNAL_NETCDF "Use internal netcdf" ON)
OPTION(INTERNAL_HDF5 "Use internal HDF5" ON)
OPTION(INTERNAL_GENGETOPT "Use internal gengetopt" ON)
OPTION(INTERNAL_QHULL "Use internal qhull" ON)

IF(INTERNAL_GENGETOPT)
  SET(GENGETOPT_URL "ftp://ftp.gnu.org/gnu/gengetopt/gengetopt-2.22.5.tar.gz" CACHE STRING "URL to download gengetopt from")
  mark_as_advanced(GENGETOPT_URL)
ENDIF(INTERNAL_GENGETOPT)

IF(INTERNAL_HDF5)
  SET(HDF5_URL "https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.10/hdf5-1.10.5/src/hdf5-1.10.5.tar.gz" CACHE STRING "URL to download HDF5 from")
  mark_as_advanced(HDF5_URL)
ENDIF(INTERNAL_HDF5)

IF(INTERNAL_NETCDF)
  SET(NETCDF_URL "https://github.com/Unidata/netcdf-c/archive/v4.7.3.tar.gz" CACHE STRING "URL to download NetCDF from")
  SET(NETCDFCXX_URL "https://github.com/Unidata/netcdf-cxx4/archive/v4.3.1.tar.gz" CACHE STRING "URL to download NetCDF-CXX from")
  mark_as_advanced(NETCDF_URL)
ENDIF(INTERNAL_NETCDF)

IF(INTERNAL_BOOST)
  SET(BOOST_URL "https://boostorg.jfrog.io/artifactory/main/release/1.76.0/source/boost_1_76_0.tar.gz" CACHE STRING "URL to download Boost from")
  mark_as_advanced(BOOST_URL)
ELSE(INTERNAL_BOOST)
  find_package(Boost 1.74.0 COMPONENTS format spirit phoenix python FATAL_ERROR)
ENDIF(INTERNAL_BOOST)

IF(INTERNAL_QHULL)
  SET(QHULL_URL "${CMAKE_SOURCE_DIR}/external/qhull-2012.1-src.tgz" CACHE STRING "URL to download QHull from")
  mark_as_advanced(QHULL_URL)
ENDIF(INTERNAL_QHULL)

SET(HEALPIX_URL "https://sourceforge.net/projects/healpix/files/Healpix_3.50/healpix_cxx-3.50.0.tar.gz/download" CACHE STRING "URL for Healpix")
mark_as_advanced(HEALPIX_URL)


find_library(ZLIB_LIBRARY z)
find_library(DL_LIBRARY dl)

SET(CONFIGURE_CPP_FLAGS "${EXTRA_CPP_FLAGS}")
SET(CONFIGURE_LD_FLAGS "$ENV{LDFLAGS}")
string(REPLACE ":" ";" OUT_LIBRARY_PATH "$ENV{LIBRARY_PATH}")
foreach(P ${OUT_LIBRARY_PATH})
  SET(CONFIGURE_LD_FLAGS "${CONFIGURE_LD_FLAGS} -L${P}")
endforeach()
message(STATUS "Extra flags for configure: ${CONFIGURE_LD_FLAGS}")


##################
# Build gengetopt
##################

if (INTERNAL_GENGETOPT)
  SET(GENGETOPT_SOURCE_DIR ${BUILD_PREFIX}/gengetopt-prefix/src/gengetopt)
  SET(GENGETOPT_BIN_DIR ${EXT_INSTALL})
  ExternalProject_Add(gengetopt
    PREFIX ${BUILD_PREFIX}/gengetopt-prefix
    URL ${GENGETOPT_URL}
    CONFIGURE_COMMAND ${GENGETOPT_SOURCE_DIR}/configure
                --prefix=${GENGETOPT_BIN_DIR}
		CPPFLAGS=${CONFIGURE_CPP_FLAGS}
		LDFLAGS=${CONFIGURE_LD_FLAGS}
		CC=${CMAKE_C_COMPILER}
		CXX=${CMAKE_CXX_COMPILER}
    INSTALL_COMMAND ${CMAKE_MAKE_PROGRAM} install
  )
  SET(GENGETOPT ${GENGETOPT_BIN_DIR}/bin/gengetopt CACHE FILEPATH "Path GenGetOpt binary")
else(INTERNAL_GENGETOPT)
  find_program(GENGETOPT gengetopt)
endif(INTERNAL_GENGETOPT)
mark_as_advanced(GENGETOPT)

###############
# Build HDF5
###############

if (INTERNAL_HDF5)

  SET(HDF5_SOURCE_DIR ${BUILD_PREFIX}/hdf5-prefix/src/hdf5)
  SET(HDF5_BIN_DIR ${EXT_INSTALL})
  ExternalProject_Add(hdf5
    PREFIX ${BUILD_PREFIX}/hdf5-prefix
    URL ${HDF5_URL}
    URL_HASH MD5=e115eeb66e944fa7814482415dd21cc4
    CMAKE_ARGS
      -DCMAKE_INSTALL_PREFIX=${EXT_INSTALL}
      -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
      -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
      -DHDF5_BUILD_CPP_LIB=ON
      -DHDF5_BUILD_TOOLS=ON
      -DHDF5_BUILD_HL_LIB=ON
      -DBUILD_SHARED_LIBS=ON
      BUILD_BYPRODUCTS ${EXT_INSTALL}/lib/libhdf5-static.a ${EXT_INSTALL}/lib/libhdf5_cpp.a
  )
  SET(hdf5_built hdf5)
  set(HDF5_LIBRARIES ${HDF5_BIN_DIR}/lib/libhdf5.a CACHE STRING "HDF5 lib" FORCE)
  set(HDF5_HL_LIBRARIES ${HDF5_BIN_DIR}/lib/libhdf5_hl.a CACHE STRING "HDF5 lib" FORCE)
  set(HDF5_CXX_LIBRARIES ${HDF5_BIN_DIR}/lib/libhdf5_cpp.a CACHE STRING "HDF5 C++ lib" FORCE)
  SET(HDF5_INCLUDE_DIR ${HDF5_BIN_DIR}/include CACHE STRING "HDF5 include path" FORCE)
  SET(CONFIGURE_LDFLAGS "${CONFIGURE_LDFLAGS} -L${HDF5_BIN_DIR}/lib")
  SET(HDF5_ROOTDIR ${HDF5_BIN_DIR})
  SET(ares_DEPS ${ares_DEPS} hdf5)
  mark_as_advanced(HDF5_LIBRARIES HDF5_HL_LIBRARIES HDF5_CXX_LIBRARIES HDF5_INCLUDE_DIR)
else(INTERNAL_HDF5)
  mark_as_advanced(CLEAR HDF5_LIBRARIES HDF5_CXX_LIBRARIES HDF5_INCLUDE_DIR)
  find_package(HDF5 COMPONENTS CXX)
  SET(HDF5_ROOTDIR ${HDF5_BIN_DIR})
  SET(HDF5_INCLUDE_DIR ${HDF5_INCLUDE_DIRS})
endif (INTERNAL_HDF5)

IF(HDF5_INCLUDE_DIR)
  SET(CONFIGURE_CPP_FLAGS "${CONFIGURE_CPP_FLAGS} -I${HDF5_INCLUDE_DIR}")
ENDIF()
SET(CONFIGURE_LIBS_FLAGS "${DL_LIB}")
mark_as_advanced(HDF5_INCLUDE_PATH HDF5_LIBRARIES HDF5_HL_LIBRARIES HDF5_CXX_LIBRARIES)

###############
# Build NetCDF
###############


if (INTERNAL_NETCDF)
  SET(NETCDF_SOURCE_DIR ${BUILD_PREFIX}/netcdf-prefix/src/netcdf)
  SET(NETCDF_CXX_SOURCE_DIR ${BUILD_PREFIX}/netcdf-cxx-prefix/src/netcdf_cxx)
  SET(NETCDF_BIN_DIR ${EXT_INSTALL})
  SET(CONFIGURE_CPP_FLAGS "${CONFIGURE_CPP_FLAGS} -I${NETCDF_BIN_DIR}/include")
  SET(CONFIGURE_LD_FLAGS "${CONFIGURE_LD_FLAGS} -L${NETCDF_BIN_DIR}/lib")
  SET(EXTRA_NC_FLAGS CPPFLAGS=${CONFIGURE_CPP_FLAGS} LDFLAGS=${CONFIGURE_LD_FLAGS} LIBS=${CONFIGURE_LIBS_FLAGS})
  ExternalProject_Add(netcdf
    DEPENDS ${hdf5_built}
    PREFIX ${BUILD_PREFIX}/netcdf-prefix
    URL ${NETCDF_URL}
    PATCH_COMMAND  ${CMAKE_COMMAND}
      -DBUILD_PREFIX=${BUILD_PREFIX}/netcdf-prefix
      -DPATCH_FILE=${CMAKE_SOURCE_DIR}/external/patch_netcdf
      -DSOURCE_PREFIX=${BUILD_PREFIX}/netcdf-prefix/src/netcdf/ncgen3
      -P ${CMAKE_SOURCE_DIR}/external/check_and_apply_patch.cmake
      CONFIGURE_COMMAND env PATH=${EXT_INSTALL}/bin:$ENV{PATH} ${NETCDF_SOURCE_DIR}/configure
         --prefix=${NETCDF_BIN_DIR} --libdir=${NETCDF_BIN_DIR}/lib
         --enable-netcdf-4  --with-pic --disable-shared --disable-dap
         --disable-cdmremote --disable-rpc
         --disable-examples ${EXTRA_NC_FLAGS} CC=${CMAKE_C_COMPILER}
         CXX=${CMAKE_CXX_COMPILER}
    BUILD_IN_SOURCE 1
    INSTALL_COMMAND ${CMAKE_MAKE_PROGRAM} install
  )

  ExternalProject_Add(netcdf_cxx
    DEPENDS netcdf
    PREFIX ${BUILD_PREFIX}/netcdf-cxx-prefix
    URL ${NETCDFCXX_URL}
    CONFIGURE_COMMAND env PATH=${EXT_INSTALL}/bin:$ENV{PATH} ${NETCDF_CXX_SOURCE_DIR}/configure
      --prefix=${NETCDF_BIN_DIR} --libdir=${NETCDF_BIN_DIR}/lib
      --with-pic --disable-shared
      --disable-examples ${EXTRA_NC_FLAGS} CC=${CMAKE_C_COMPILER}
      CXX=${CMAKE_CXX_COMPILER}
    BUILD_IN_SOURCE 1
    INSTALL_COMMAND ${CMAKE_MAKE_PROGRAM} install
  )

  SET(CONFIGURE_CPP_LDFLAGS "${CONFIGURE_LDFLAGS}")
  SET(EXTRA_NC_FLAGS CPPFLAGS=${CONFIGURE_CPP_FLAGS} LDFLAGS=${CONFIGURE_CPP_LDFLAGS})
  SET(cosmotool_DEPS ${cosmotool_DEPS} netcdf netcdf_cxx)
  SET(NETCDF_LIBRARY ${NETCDF_BIN_DIR}/lib/libnetcdf.a CACHE STRING "NetCDF lib" FORCE)
  SET(NETCDFCPP_LIBRARY ${NETCDF_BIN_DIR}/lib/libnetcdf_c++4.a CACHE STRING "NetCDF-C++ lib" FORCE)
  SET(NETCDF_INCLUDE_PATH ${NETCDF_BIN_DIR}/include CACHE STRING "NetCDF include" FORCE)
  SET(NETCDFCPP_INCLUDE_PATH ${NETCDF_INCLUDE_PATH} CACHE STRING "NetCDF C++ include path" FORCE)

ELSE(INTERNAL_NETCDF)
  find_library(NETCDF_LIBRARY netcdf)
  find_library(NETCDFCPP_LIBRARY netcdf_c++)
  find_path(NETCDF_INCLUDE_PATH NAMES netcdf.h)
  find_path(NETCDFCPP_INCLUDE_PATH NAMES netcdf)
  SET(CONFIGURE_CPP_FLAGS ${CONFIGURE_CPP_FLAGS}
          -I${NETCDF_INCLUDE_PATH} -I${NETCDFCPP_INCLUDE_PATH})
endif (INTERNAL_NETCDF)
mark_as_advanced(NETCDF_LIBRARY NETCDFCPP_LIBRARY NETCDF_INCLUDE_PATH NETCDFCPP_INCLUDE_PATH)

##################
# Build BOOST
##################

if (INTERNAL_BOOST)
  SET(cosmotool_DEPS ${cosmotool_DEPS} boost)
  SET(BOOST_SOURCE_DIR ${BUILD_PREFIX}/boost-prefix/src/boost)

  set(LINKER_EXTRA_FLAGS)
  message(STATUS "Compiler version is ${CMAKE_CXX_COMPILER_VERSION}")
  string(REGEX REPLACE "^([0-9]+\\.[0-9]+).*$" "\\1" ToolsetVer "${CMAKE_CXX_COMPILER_VERSION}")
  IF(CMAKE_CXX_COMPILER_ID MATCHES "^Intel$")
     SET(b2_toolset intel)
     SET(COMPILER_EXTRA_FLAGS "-fPIC")
  elseif (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
     if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
       SET(b2_toolset darwin)
     else()
       SET(b2_toolset gcc)
       SET(COMPILER_EXTRA_FLAGS "-fPIC -std=gnu++14")
     endif()
     add_definitions("-Wno-unused-local-typedefs")
  elseif (CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
     SET(b2_toolset darwin)
  elseif (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
     SET(b2_toolset clang)
  endif()
  SET(COMPILER_EXTRA_FLAGS "${COMPILER_EXTRA_FLAGS} -I${EXT_INSTALL}/include")
  SET(LINKER_EXTRA_FLAGS "${LINKER_EXTRA_FLAGS} -L${EXT_INSTALL}/lib")

  message(STATUS "Building boost with toolset ${b2_toolset}")

  SET(BOOST_LIBRARIES
    ${EXT_INSTALL}/lib/libboost_exception.a
  )

  ExternalProject_Add(boost
    URL ${BOOST_URL}
    URL_HASH SHA1=a5ab6eaf31d1ca181a17ecffef9d58d40d87c71d
    PREFIX ${BUILD_PREFIX}/boost-prefix
    CONFIGURE_COMMAND ${CMAKE_COMMAND}  -DTOOLSET=${b2_toolset} "-DCOMPILER:STRING=${CMAKE_CXX_COMPILER}" "-DCOMPILER_EXTRA_FLAGS=${COMPILER_EXTRA_FLAGS}" "-DINSTALL_PATH:STRING=${EXT_INSTALL}" "-DLINKER_EXTRA_FLAGS=${LINKER_EXTRA_FLAGS}" "-DSRC_DIR:STRING=${BOOST_SOURCE_DIR}" -P ${CMAKE_SOURCE_DIR}/external/configure_boost.cmake
    BUILD_IN_SOURCE 1
    BUILD_COMMAND ${BOOST_SOURCE_DIR}/b2 --with-exception toolset=${b2_toolset}-cmake variant=release
    INSTALL_COMMAND  ${BOOST_SOURCE_DIR}/b2 install --with-exception
    BUILD_BYPRODUCTS ${BOOST_LIBRARIES}
  )
  set(Boost_INCLUDE_DIRS ${BOOST_SOURCE_DIR} CACHE STRING "Boost path" FORCE)
  set(Boost_INTERNAL_INSTALL ${EXT_INSTALL})
  ##  set(Boost_LIBRARIES ${BOOST_SOURCE_DIR}/stage/lib/libboost_python.a CACHE STRING "Boost libraries" FORCE)
endif (INTERNAL_BOOST)
mark_as_advanced(Boost_INCLUDE_DIRS Boost_LIBRARIES)

##################
# Build GSl
##################

IF(INTERNAL_GSL)
  SET(GSL_URL "http://ftpmirror.gnu.org/gsl/gsl-2.3.tar.gz" CACHE STRING "URL to download GSL from ")
  mark_as_advanced(GSL_URL)

  SET(GSL_SOURCE_DIR ${BUILD_PREFIX}/gsl-prefix/src/gsl)
  ExternalProject_Add(gsl
    URL ${GSL_URL}
    PREFIX ${BUILD_PREFIX}/gsl-prefix
    CONFIGURE_COMMAND ${GSL_SOURCE_DIR}/configure
           --prefix=${EXT_INSTALL} --disable-shared
           --with-pic --libdir=${EXT_INSTALL}/lib
           CPPFLAGS=${CONFIGURE_CPP_FLAGS} CC=${CMAKE_C_COMPILER} CXX=${CMAKE_CXX_COMPILER}
    BUILD_IN_SOURCE 1
    BUILD_BYPRODUCTS ${EXT_INSTALL}/lib/libgsl.a ${EXT_INSTALL}/lib/libgslcblas.a
  )
  SET(GSL_INTERNAL_LIBS ${EXT_INSTALL}/lib)
  SET(GSL_LIBRARY ${GSL_INTERNAL_LIBS}/libgsl.a CACHE STRING "GSL internal path" FORCE)
  SET(GSL_CBLAS_LIBRARY ${GSL_INTERNAL_LIBS}/libgslcblas.a CACHE STRING "GSL internal path" FORCE)
  set(GSL_INCLUDE ${CMAKE_BINARY_DIR}/ext_install/include CACHE STRING "GSL internal path" FORCE)
  SET(cosmotool_DEPS ${cosmotool_DEPS} gsl)
  mark_as_advanced(GSL_LIBRARY GSL_INCLUDE GSL_CBLAS_LIBRARY)
ELSE(INTERNAL_GSL)
  mark_as_advanced(CLEAR GSL_LIBRARY GSL_INCLUDE GSL_CBLAS_LIBRARY)
  if (NOT DEFINED GSL_LIBRARY OR NOT GSL_LIBRARY)
    find_library(GSL_LIBRARY gsl)
  endif()
  if (NOT DEFINED GSL_CBLAS_LIBRARY OR NOT GSL_CBLAS_LIBRARY)
    find_library(GSL_CBLAS_LIBRARY gslcblas)
  endif()
  if (NOT DEFINED GSL_INCLUDES OR NOT GSL_INCLUDES)
    find_path(GSL_INCLUDE NAMES gsl/gsl_blas.h )
  endif()
  if (NOT (GSL_LIBRARY OR GSL_CBLAS_LIBRARY OR GSL_INCLUDE))
    message(FATAL_ERROR "GSL has not been found")
  endif()
ENDIF(INTERNAL_GSL)
message(STATUS "GSL paths: ${GSL_LIBRARY} ${GSL_CBLAS_LIBRARY} ${GSL_INCLUDE}")


#################
# Build cfitsio
#################

ExternalProject_Add(cfitsio
  URL file://${CMAKE_SOURCE_DIR}/external/cfitsio-3.47.tar.gz
  URL_HASH SHA1=5a25016dcaf12117d950e4278e10d39c6c7d33a5
  PREFIX ${BUILD_PREFIX}/cfitsio-prefix
  CONFIGURE_COMMAND ./configure --prefix=${EXT_INSTALL} --disable-curl --libdir=${EXT_INSTALL}/lib CPPFLAGS=${CONFIGURE_CPP_FLAGS} CC=${CMAKE_C_COMPILER} CXX=${CMAKE_CXX_COMPILER}
  BUILD_IN_SOURCE 1
  BUILD_BYPRODUCTS ${EXT_INSTALL}/lib/libcfitsio.a
)

SET(CFITSIO_LIBRARY ${EXT_INSTALL}/lib/libcfitsio.a)
SET(CFITSIO_INCLUDE_PATH ${EXT_INSTALL}/include)

#################
# Build Healpix
#################
SET(HEALPIX_BUILD ${BUILD_PREFIX}/healpix-prefix/src/healpix-build)
SET(HEALPIX_DIR ${BUILD_PREFIX}/healpix-prefix/src/healpix)

ExternalProject_Add(healpix
    DEPENDS cfitsio
    PREFIX ${BUILD_PREFIX}/healpix-prefix
    URL ${HEALPIX_URL}
    URL_HASH SHA1=c8a537e743f760dfa453cad246065d37f72fc0cb
    CONFIGURE_COMMAND ${CMAKE_COMMAND}
        -DHEALPIX_CC=${CMAKE_C_COMPILER}
        -DHEALPIX_CXX=${CMAKE_CXX_COMPILER}
        -DHEALPIX_DIR:STRING=${HEALPIX_DIR}
        -DHEALPIX_INSTALL:STRING=${EXT_INSTALL}
        -DCFITSIO_LIB:STRING=${CFITSIO_LIBRARY}
        -P ${CMAKE_SOURCE_DIR}/external/configure_healpix.cmake
)
set(HPIX_LIBPATH ${EXT_INSTALL}/lib)
set(HEALPIX_LIBRARY ${HPIX_LIBPATH}/libhealpix_cxx.a)

SET(HEALPIX_INCLUDE_PATH ${EXT_INSTALL}/include/healpix_cxx)
SET(HEALPIX_LIBRARIES ${HEALPIX_LIBRARY} ${CFITSIO_LIBRARY} )
set(GSL_LIBRARIES ${GSL_LIBRARY} ${GSL_CBLAS_LIBRARY})
SET(NETCDF_LIBRARIES ${NETCDFCPP_LIBRARY} ${NETCDF_LIBRARY} ${HDF5_HL_LIBRARIES} ${HDF5_LIBRARIES} ${ZLIB_LIBRARY})

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
  message(FATAL_ERROR "Only packaged QHull is supported")
endif(INTERNAL_QHULL)

SET(QHULL_LIBRARIES ${QHULL_CPP_LIBRARY} ${QHULL_LIBRARY} )


###############
# Build libSDF
###############
IF(SDF_SUPPORT)
  SET(LIBSDF_ARCH x86_64 CACHE STRING "SDF architecture to activate")
  mark_as_advanced(LIBSDF_ARCH)
  SET(LIBSDF_PATH ${CMAKE_SOURCE_DIR}/external/libsdf)

  ExternalProject_Add(libSDF
    URL ${CMAKE_SOURCE_DIR}/external/mswarren-libsdf-b4b9f9464b5b.tar.gz
    PREFIX ${BUILD_PREFIX}/libSDF-prefix
    CONFIGURE_COMMAND echo No configure
    BUILD_COMMAND ${CMAKE_MAKE_PROGRAM} -f Make.simple
    INSTALL_COMMAND ${CMAKE_COMMAND} -DDEST_DIR=${CMAKE_BINARY_DIR}/ext_build/sdf -DLIBSDF_ARCH=${LIBSDF_ARCH} -DLIBSDF_PATH=${BUILD_PREFIX}/libSDF-prefix/src/libSDF -P ${CMAKE_SOURCE_DIR}/external/install_sdf.cmake
    BUILD_IN_SOURCE 1
    PATCH_COMMAND  ${CMAKE_COMMAND}
      -DBUILD_PREFIX=${BUILD_PREFIX}/libSDF-prefix
      -DPATCH_FILE=${CMAKE_SOURCE_DIR}/external/patch_sdf
      -DSOURCE_PREFIX=${BUILD_PREFIX}/libSDF-prefix/src/libSDF
      -P ${CMAKE_SOURCE_DIR}/external/check_and_apply_patch.cmake
  )
  SET(LIBSDF_INCLUDE_PATH ${BUILD_PREFIX}/libSDF-prefix/src)
  SET(LIBSDF_LIBRARY ${BUILD_PREFIX}/libSDF-prefix/src/libSDF/libSDF.a)

  find_library(RT_LIBRARY rt)
  IF (RT_LIBRARY)
    SET(LIBSDF_LIBRARY ${LIBSDF_LIBRARY} ${RT_LIBRARY})
  ENDIF (RT_LIBRARY)
ENDIF(SDF_SUPPORT)

include_directories(${CMAKE_BINARY_DIR}/src
                    ${NETCDF_INCLUDE_PATH} ${GSL_INCLUDE_PATH}
                    ${HDF5_INCLUDE_PATH}
                    ${Boost_INCLUDE_DIRS}
		    ${QHULL_INCLUDE_PATH} ${LIBSDF_INCLUDE_PATH})

message(STATUS "Boost_INTERNAL_INSTALL=${Boost_INTERNAL_INSTALL}")
