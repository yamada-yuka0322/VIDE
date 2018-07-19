include(FindOpenMP)

OPTION(ENABLE_OPENMP "Set to Yes if Healpix and/or you need openMP" OFF)

SET(FFTW_URL "http://www.fftw.org/fftw-3.3.3.tar.gz" CACHE URL "URL to download FFTW from")
SET(EIGEN_URL "http://bitbucket.org/eigen/eigen/get/3.2.10.tar.gz" CACHE URL "URL to download Eigen from")
SET(GENGETOPT_URL "ftp://ftp.gnu.org/gnu/gengetopt/gengetopt-2.22.5.tar.gz" CACHE STRING "URL to download gengetopt from")
SET(HDF5_URL "https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.8/hdf5-1.8.18/src/hdf5-1.8.18.tar.bz2" CACHE STRING "URL to download HDF5 from")
SET(NETCDF_URL "ftp://ftp.unidata.ucar.edu/pub/netcdf/netcdf-4.5.0.tar.gz" CACHE STRING "URL to download NetCDF from")
SET(NETCDFCXX_URL "https://github.com/Unidata/netcdf-cxx4/archive/v4.3.0.tar.gz" CACHE STRING "URL to download NetCDF-C++ from")
SET(BOOST_URL "http://sourceforge.net/projects/boost/files/boost/1.61.0/boost_1_61_0.tar.gz/download" CACHE STRING "URL to download Boost from")
SET(GSL_URL "ftp://ftp.gnu.org/gnu/gsl/gsl-1.15.tar.gz" CACHE STRING "URL to download GSL from ")
mark_as_advanced(FFTW_URL EIGEN_URL HDF5_URL NETCDF_URL BOOST_URL GSL_URL)


MACRO(CHECK_CHANGE_STATE VAR) 
  IF (DEFINED _PREVIOUS_${VAR})
    IF (NOT ${_PREVIOUS_${VAR}} EQUAL ${${VAR}})
      foreach(loopvar ${ARGN})
         UNSET(${loopvar} CACHE)
      endforeach()
    ENDIF (NOT ${_PREVIOUS_${VAR}} EQUAL ${${VAR}})
  ENDIF (DEFINED _PREVIOUS_${VAR})
  SET(_PREVIOUS_${VAR} ${${VAR}} CACHE INTERNAL "Internal value")
ENDMACRO(CHECK_CHANGE_STATE)

CHECK_CHANGE_STATE(INTERNAL_BOOST Boost_LIBRARIES Boost_INCLUDE_DIRS)
CHECK_CHANGE_STATE(INTERNAL_EIGEN EIGEN3_INCLUDE_DIRS)
CHECK_CHANGE_STATE(INTERNAL_GSL GSL_LIBRARY GSL_CBLAS_LIBRARY GSL_INCLUDE)
CHECK_CHANGE_STATE(INTERNAL_HDF5 
    HDF5_INCLUDE_DIR HDF5_LIBRARIES HDF5_CXX_LIBRARIES 
    HDF5_C_STATIC_LIBRARY HDF5_HL_STATIC_LIBRARY HDF5_CXX_STATIC_LIBRARY)
CHECK_CHANGE_STATE(INTERNAL_DLIB DLIB_INCLUDE_DIR DLIB_LIBRARIES)


IF(ENABLE_OPENMP)
  IF (NOT OPENMP_FOUND)
    MESSAGE(ERROR "No known compiler option for enabling OpenMP")
  ENDIF(NOT OPENMP_FOUND)

  SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_C_FLAGS}")
ENDIF(ENABLE_OPENMP)


SET(BUILD_PREFIX ${CMAKE_BINARY_DIR}/external_build)
SET(EXT_INSTALL ${CMAKE_BINARY_DIR}/ext_install)
SET(CONFIGURE_LIBS )
SET(CONFIGURE_CPP_FLAGS "")
SET(CONFIGURE_LDFLAGS "")

if (ENABLE_SHARP)
  SET(DEP_BUILD ${BUILD_PREFIX}/sharp-prefix/src/sharp/auto)
  IF(NOT ENABLE_OPENMP)
    SET(SHARP_OPENMP --disable-openmp)
  ENDIF()
  ExternalProject_Add(sharp
    URL ${CMAKE_SOURCE_DIR}/external/libsharp-6077806.tar.gz
    PREFIX ${BUILD_PREFIX}/sharp-prefix
    BUILD_IN_SOURCE 1 
    CONFIGURE_COMMAND autoconf && ./configure "CC=${CMAKE_C_COMPILER}" "CXX=${CMAKE_CXX_COMPILER}" --prefix=${DEP_BUILD} ${SHARP_OPENMP}
    BUILD_COMMAND ${CMAKE_MAKE_PROGRAM}
    INSTALL_COMMAND echo "No install"
  )
  SET(CUTILS_LIBRARY ${DEP_BUILD}/lib/libc_utils.a)
  SET(FFTPACK_LIBRARY ${DEP_BUILD}/lib/libfftpack.a)
  SET(SHARP_LIBRARY ${DEP_BUILD}/lib/libsharp.a)
  SET(SHARP_LIBRARIES ${SHARP_LIBRARY} ${FFTPACK_LIBRARY} ${CUTILS_LIBRARY})
  SET(SHARP_INCLUDE_PATH ${DEP_BUILD}/include)
endif (ENABLE_SHARP)


###############
# Build HDF5
###############

if (INTERNAL_HDF5)
  SET(HDF5_SOURCE_DIR ${BUILD_PREFIX}/hdf5-prefix/src/hdf5)
  SET(HDF5_BIN_DIR ${EXT_INSTALL})
  ExternalProject_Add(hdf5
    PREFIX ${BUILD_PREFIX}/hdf5-prefix
    URL ${HDF5_URL}
    URL_HASH MD5=29117bf488887f89888f9304c8ebea0b
    CMAKE_ARGS
       -DCMAKE_INSTALL_PREFIX=${EXT_INSTALL}
       -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
       -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
       -DHDF5_BUILD_CPP_LIB=ON
       -DHDF5_BUILD_TOOLS=ON
       -DHDF5_BUILD_HL_LIB=ON
       -DBUILD_SHARED_LIBS=OFF
  )
  SET(cosmotool_DEPS ${cosmotool_DEPS} hdf5)
  SET(hdf5_built hdf5)
  SET(ENV{HDF5_ROOT} ${HDF5_BIN_DIR})
  SET(HDF5_ROOTDIR ${HDF5_BIN_DIR})
  SET(CONFIGURE_LDFLAGS "${CONFIGURE_LDFLAGS} -L${HDF5_BIN_DIR}/lib")
  SET(CONFIGURE_LIBS "${CONFIGURE_LIBS} -ldl")
  set(HDF5_C_STATIC_LIBRARY ${HDF5_BIN_DIR}/lib/libhdf5-static.a)
  set(HDF5_HL_STATIC_LIBRARY ${HDF5_BIN_DIR}/lib/libhdf5_hl-static.a)
  set(HDF5_LIBRARIES ${HDF5_BIN_DIR}/lib/libhdf5-static.a CACHE STRING "HDF5 lib" FORCE)
  set(HDF5_HL_LIBRARIES ${HDF5_BIN_DIR}/lib/libhdf5_hl-static.a CACHE STRING "HDF5 HL lib" FORCE)
  set(HDF5_CXX_LIBRARIES ${HDF5_BIN_DIR}/lib/libhdf5_cpp-static.a CACHE STRING "HDF5 C++ lib" FORCE)
  SET(HDF5_INCLUDE_DIRS ${HDF5_BIN_DIR}/include CACHE STRING "HDF5 include path" FORCE)
  mark_as_advanced(HDF5_LIBRARIES HDF5_CXX_LIBRARIES HDF5_INCLUDE_DIRS)

  MESSAGE(STATUS "Internal HDF5 directory: $ENV{HDF5_ROOT}")
  MESSAGE(STATUS "Libs: ${HDF5_LIBRARIES}")
  SET(HDF5_FOUND TRUE)
else (INTERNAL_HDF5)
  mark_as_advanced(CLEAR HDF5_LIBRARIES HDF5_CXX_LIBRARIES HDF5_INCLUDE_DIRS)
  if(HDF5_ROOTDIR)
    SET(ENV{HDF5_ROOT} ${HDF5_ROOTDIR})
  endif(HDF5_ROOTDIR)
  find_package(HDF5 CONFIG QUIET COMPONENTS C CXX HL static)
  if (NOT HDF5_FOUND)
    cmessage(CWARNING "Could not find HDF5 cmake config. Try classical exploration")
    find_package(HDF5 COMPONENTS C CXX HL)
    cmessage(STATUS "HDF5 lib: ${HDF5_LIBRARIES}")
    cmessage(STATUS "HDF5 includes: ${HDF5_INCLUDE_DIRS}")
    cmessage(STATUS "HDF5 C lib: ${HDF5_C_LIBRARY}")
    cmessage(STATUS "HDF5 HL lib: ${HDF5_HL_LIBRARY}")
    cmessage(STATUS "HDF5 BIN: ${HDF5_BIN_DIR}")
    foreach(hdf5lib IN LISTS HDF5_LIBRARIES)
      if (${hdf5lib} MATCHES "(hdf5)|(HDF5)")
        get_filename_component(HDF5_BIN_DIR ${hdf5lib} DIRECTORY)
      endif()
    endforeach()
    cmessage(STATUS "HDF5 libpath: ${HDF5_BIN_DIR}")
  else()
    cmessage(STATUS "Found HDF5 cmake config.")
    cmessage(STATUS "HDF5_C_STATIC_LIBRARY : ${HDF5_C_STATIC_LIBRARY}")
    set(HDF5_LIBRARIES ${HDF5_C_STATIC_LIBRARY} CACHE STRING "HDF5 lib" FORCE)
    set(HDF5_HL_LIBRARIES ${HDF5_HL_STATIC_LIBRARY} CACHE STRING "HDF5 HL lib" FORCE)
    set(HDF5_CXX_LIBRARIES ${HDF5_CXX_STATIC_LIBRARY} CACHE STRING "HDF5 C++ lib" FORCE)
    get_filename_component(HDF5_BIN_DIR ${HDF5_C_STATIC_LIBRARY} DIRECTORY)
  endif()
  SET(CONFIGURE_LDFLAGS "${CONFIGURE_LDFLAGS} -L${HDF5_BIN_DIR}")
endif (INTERNAL_HDF5)

foreach(include_dir ${HDF5_INCLUDE_DIRS})
  SET(CONFIGURE_CPP_FLAGS "${CONFIGURE_CPP_FLAGS} -I${include_dir}")
endforeach(include_dir)

###############
# Build NetCDF
###############


if (INTERNAL_NETCDF)
  SET(NETCDF_SOURCE_DIR ${BUILD_PREFIX}/netcdf-prefix/src/netcdf)
  SET(NETCDF_BIN_DIR ${EXT_INSTALL})
  SET(CONFIGURE_CPP_FLAGS "${CONFIGURE_CPP_FLAGS} -I${NETCDF_BIN_DIR}/include")
  SET(CONFIGURE_LDFLAGS "${CONFIGURE_LDFLAGS} -L${NETCDF_BIN_DIR}/lib")
  SET(EXTRA_NC_FLAGS CPPFLAGS=${CONFIGURE_CPP_FLAGS} LIBS=${CONFIGURE_LIBS} LDFLAGS=${CONFIGURE_LDFLAGS})
  SET(NETCDF_CONFIG_COMMAND ${NETCDF_SOURCE_DIR}/configure
         --prefix=${NETCDF_BIN_DIR} --libdir=${NETCDF_BIN_DIR}/lib
         --enable-netcdf-4  --with-pic --disable-shared --disable-dap 
         --disable-cdmremote --disable-rpc --enable-cxx-4 
         --disable-examples ${EXTRA_NC_FLAGS} CC=${CMAKE_C_COMPILER}
         CXX=${CMAKE_CXX_COMPILER})
  list(INSERT CMAKE_PREFIX_PATH 0 ${EXT_INSTALL})
  string(REPLACE ";" "|" CMAKE_PREFIX_PATH_ALT_SEP "${CMAKE_PREFIX_PATH}")
  ExternalProject_Add(netcdf
    DEPENDS ${hdf5_built}
    PREFIX ${BUILD_PREFIX}/netcdf-prefix
    URL ${NETCDF_URL}
  LIST_SEPARATOR |
	CMAKE_ARGS
		-DCMAKE_PREFIX_PATH=${CMAKE_PREFIX_PATH_ALT_SEP}
		-DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
		-DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
		-DBUILD_SHARED_LIBS=OFF
		-DBUILD_TESTING=OFF
		-DCMAKE_BUILD_TYPE=Release
		-DENABLE_NETCDF4=ON
    -DCMAKE_POSITION_INDEPENDENT_CODE=ON
    -DENABLE_DAP=OFF
		-DCMAKE_INSTALL_PREFIX=${NETCDF_BIN_DIR}
                -DHDF5_C_LIBRARY=${HDF5_C_STATIC_LIBRARY}
                -DHDF5_HL_LIBRARY=${HDF5_HL_STATIC_LIBRARY}
                -DHDF5_INCLUDE_DIR=${HDF5_INCLUDE_DIRS}
    -DCMAKE_INSTALL_LIBDIR=lib
  )
  
  SET(NETCDFCXX_SOURCE_DIR ${BUILD_PREFIX}/netcdf-c++-prefix/src/netcdf-c++)
  ExternalProject_Add(netcdf-c++
    DEPENDS ${hdf5_built} netcdf
    PREFIX ${BUILD_PREFIX}/netcdf-c++-prefix
    URL ${NETCDFCXX_URL}
	CMAKE_ARGS
		-DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
		-DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
		-DBUILD_SHARED_LIBS=OFF
    -DCMAKE_POSITION_INDEPENDENT_CODE=ON
		-DBUILD_TESTING=OFF
		-DCMAKE_BUILD_TYPE=Release
		-DCMAKE_INSTALL_PREFIX=${NETCDF_BIN_DIR}
    -DCMAKE_INSTALL_LIBDIR=lib
  )
#  SET(CONFIGURE_CPP_LDFLAGS "${CONFIGURE_LDFLAGS}")
#  SET(EXTRA_NC_FLAGS CPPFLAGS=${CONFIGURE_CPP_FLAGS} LDFLAGS=${CONFIGURE_CPP_LDFLAGS})
  SET(cosmotool_DEPS ${cosmotool_DEPS} netcdf netcdf-c++)
#  find_library(NETCDF_LIBRARY netcdf NO_DEFAULT_PATH HINTS ${NETCDF_BIN_DIR}/lib ${NETCDF_BIN_DIR}/lib32 ${NETCDF_BIN_DIR}/lib64)
#  find_library(NETCDFCPP_LIBRARY netcdf-cxx4 NO_DEFAULT_PATH HINTS ${NETCDF_BIN_DIR}/lib ${NETCDF_BIN_DIR}/lib32 ${NETCDF_BIN_DIR}/lib64)
  SET(NETCDF_LIBRARY ${NETCDF_BIN_DIR}/lib/libnetcdf.a CACHE STRING "NetCDF lib" FORCE)
  SET(NETCDFCPP_LIBRARY ${NETCDF_BIN_DIR}/lib/libnetcdf-cxx4.a CACHE STRING "NetCDF-C++ lib" FORCE)
  SET(NETCDF_INCLUDE_PATH ${NETCDF_BIN_DIR}/include CACHE STRING "NetCDF include" FORCE)
  SET(NETCDFCPP_INCLUDE_PATH ${NETCDF_INCLUDE_PATH} CACHE STRING "NetCDF C++ include path" FORCE)

ELSE(INTERNAL_NETCDF)
  find_path(NETCDF_INCLUDE_PATH NAMES netcdf.h)
  find_path(NETCDFCPP_INCLUDE_PATH NAMES netcdfcpp.h netcdf)
  find_library(NETCDF_LIBRARY netcdf)
  find_library(NETCDFCPP_LIBRARY NAMES netcdf_c++4 netcdf_c++)

  SET(CONFIGURE_CPP_FLAGS "${CONFIGURE_CPP_FLAGS} -I${NETCDF_INCLUDE_PATH} -I${NETCDFCPP_INCLUDE_PATH}")          
endif (INTERNAL_NETCDF)
mark_as_advanced(NETCDF_LIBRARY NETCDFCPP_LIBRARY NETCDF_INCLUDE_PATH NETCDFCPP_INCLUDE_PATH)

##################
# Build BOOST
##################

if (INTERNAL_BOOST)
  message(STATUS "Building Boost")
  SET(BOOST_SOURCE_DIR ${BUILD_PREFIX}/boost-prefix/src/boost)
  ExternalProject_Add(boost
    URL ${BOOST_URL}
    PREFIX ${BUILD_PREFIX}/boost-prefix
    CONFIGURE_COMMAND 
           ${BOOST_SOURCE_DIR}/bootstrap.sh --prefix=${CMAKE_BINARY_DIR}/ext_build/boost
    BUILD_IN_SOURCE 1
    BUILD_COMMAND ${BOOST_SOURCE_DIR}/b2 --with-exception 
    INSTALL_COMMAND echo "No install"
  )
  set(Boost_INCLUDE_DIRS ${BOOST_SOURCE_DIR} CACHE STRING "Boost path" FORCE)
  set(Boost_LIBRARIES ${BOOST_SOURCE_DIR}/stage/lib/libboost_python.a CACHE STRING "Boost libraries" FORCE)
  set(Boost_FOUND YES)
  set(Boost_DEP boost)

ELSE (INTERNAL_BOOST)
  find_package(Boost 1.53 QUIET)
  set(Boost_DEP)
  if (NOT Boost_FOUND)
    cmessage(CWARNING "Boost >= 1.53 was not found")
  endif()
endif (INTERNAL_BOOST)
mark_as_advanced(Boost_INCLUDE_DIRS Boost_LIBRARIES)

##################
# Build GSL
##################

IF(INTERNAL_GSL)
  SET(GSL_SOURCE_DIR ${BUILD_PREFIX}/gsl-prefix/src/gsl)
  ExternalProject_Add(gsl
    URL ${GSL_URL}
    PREFIX ${BUILD_PREFIX}/gsl-prefix
    CONFIGURE_COMMAND ${GSL_SOURCE_DIR}/configure
           --prefix=${EXT_INSTALL} --disable-shared
           --with-pic
           CPPFLAGS=${CONFIGURE_CPP_FLAGS} CC=${CMAKE_C_COMPILER} CXX=${CMAKE_CXX_COMPILER}
    BUILD_IN_SOURCE 1
    BUILD_COMMAND ${CMAKE_MAKE_PROGRAM}
    INSTALL_COMMAND ${CMAKE_MAKE_PROGRAM} install
  )
  SET(GSL_INTERNAL_LIBS ${EXT_INSTALL}/lib)
  SET(GSL_LIBRARY ${GSL_INTERNAL_LIBS}/libgsl.a CACHE STRING "GSL internal path" FORCE)
  SET(GSLCBLAS_LIBRARY ${GSL_INTERNAL_LIBS}/libgslcblas.a CACHE STRING "GSL internal path" FORCE)
  set(GSL_INCLUDE_PATH ${CMAKE_BINARY_DIR}/ext_build/gsl/include CACHE STRING "GSL internal path" FORCE)
  set(GSL_LIBRARIES ${GSL_LIBRARY} ${GSLCBLAS_LIBRARY})
  SET(cosmotool_DEPS ${cosmotool_DEPS} gsl)
ELSE(INTERNAL_GSL)
  find_path(GSL_INCLUDE_PATH NAMES gsl/gsl_blas.h)
  find_library(GSL_LIBRARY gsl)
  find_library(GSLCBLAS_LIBRARY gslcblas)

  set(GSL_LIBRARIES ${GSL_LIBRARY} ${GSLCBLAS_LIBRARY})

ENDIF(INTERNAL_GSL)
mark_as_advanced(GSL_LIBRARY GSLCBLAS_LIBRARY GSL_INCLUDE_PATH)


#############
# Build FFTW
#############

IF(INTERNAL_FFTW)
	SET(EXTRA_FFTW_CONF)
	IF(HAVE_SSE)
		SET(EXTRA_FFTW_CONF ${EXTRA_FFTW_CONF} --enable-sse)
	ENDIF(HAVE_SSE)
	IF(HAVE_SSE2)
		SET(EXTRA_FFTW_CONF ${EXTRA_FFTW_CONF} --enable-sse2)
	ENDIF(HAVE_SSE2)
	IF(HAVE_AVX)
		SET(EXTRA_FFTW_CONF ${EXTRA_FFTW_CONF} --enable-avx)
	ENDIF(HAVE_AVX)

  SET(cosmotool_DEPS ${cosmotool_DEPS} fftw)
  SET(FFTW_SOURCE ${BUILD_PREFIX}/fftw-prefix/src/fftw)
  ExternalProject_Add(fftw
     URL ${FFTW_URL}
     PREFIX ${BUILD_PREFIX}/fftw-prefix
     CONFIGURE_COMMAND 
           ${FFTW_SOURCE}/configure
                 --prefix=${EXT_INSTALL}
		 ${EXTRA_FFTW_CONF} --disable-shared --enable-threads
     BUILD_COMMAND ${CMAKE_MAKE_PROGRAM}
     INSTALL_COMMAND ${CMAKE_MAKE_PROGRAM} install
  )
  SET(FFTW3_LIBRARY_DIRS ${EXT_INSTALL}/lib)
  SET(FFTW3_INCLUDE_PATH ${EXT_INSTALL}/include)
  SET(FFTW3_THREADS ${EXT_INSTALL}/lib/libfftw3_threads.a)
  SET(FFTW3_LIBRARIES ${EXT_INSTALL}/lib/libfftw3.a)
  
ELSE (INTERNAL_FFTW)  
  pkg_check_modules(FFTW3 fftw3>=3.3)
  pkg_check_modules(FFTW3F fftw3f>=3.3)

  find_library(FFTW3F_LIBRARY_FULL fftw3f PATHS ${FFTW3F_LIBDIR} NO_DEFAULT_PATH)
  find_library(FFTW3_LIBRARY_FULL fftw3 PATHS ${FFTW3_LIBDIR} NO_DEFAULT_PATH)

ENDIF(INTERNAL_FFTW)                   


##############
# Build Eigen
##############
IF (INTERNAL_EIGEN)
  ExternalProject_Add(eigen
     URL ${EIGEN_URL}
     URL_HASH SHA256=04f8a4fa4afedaae721c1a1c756afeea20d3cdef0ce3293982cf1c518f178502 
     PREFIX ${BUILD_PREFIX}/eigen-prefix
     CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${EXT_INSTALL}
       -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
       -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
  )
  SET(EIGEN3_INCLUDE_DIRS ${EXT_INSTALL}/include/eigen3)

ELSE (INTERNAL_EIGEN)
   if(DEFINED EIGEN_PATH)
     set(_eigen_old_pkg_path $ENV{PKG_CONFIG_PATH})
     set(ENV{PKG_CONFIG_PATH} ${EIGEN_PATH}/share/pkgconfig)
   endif()
   pkg_check_modules(EIGEN3 NO_CMAKE_PATH NO_CMAKE_ENVIRONMENT_PATH REQUIRED eigen3)
   if(DEFINED EIGEN_PATH)
     set(ENV{PKG_CONFIG_PATH} ${_eigen_old_pkg_path})
   endif()
  if (NOT  EIGEN3_FOUND)
    cmessage(CWARNING "Eigen library not found")
  else()
    cmessage(STATUS "Found EIGEN3 in ${EIGEN3_INCLUDE_DIRS}")
  endif()
ENDIF(INTERNAL_EIGEN)



SET(cosmotool_DEPS ${cosmotool_DEPS} omptl)
SET(OMPTL_BUILD_DIR ${BUILD_PREFIX}/omptl-prefix/src/omptl)
ExternalProject_Add(omptl
  PREFIX ${BUILD_PREFIX}/omptl-prefix
  URL ${CMAKE_SOURCE_DIR}/external/omptl-20120422.tar.bz2
  CONFIGURE_COMMAND echo "No configure"
  BUILD_COMMAND echo "No build"
  PATCH_COMMAND patch -p1 -t -N < ${CMAKE_SOURCE_DIR}/external/patch-omptl
  INSTALL_COMMAND ${CMAKE_COMMAND} -E copy_directory ${OMPTL_BUILD_DIR} ${EXT_INSTALL}/include/omptl
)
include_directories(${EXT_INSTALL}/include)
##include_directories(${OMPTL_BUILD_DIR}/src/)

