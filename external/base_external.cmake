include(FindOpenMP)

OPTION(ENABLE_OPENMP "Set to Yes if Healpix and/or you need openMP" OFF)
OPTION(SDF_SUPPORT "Set to Yes to activate support for SDF" ON)

IF(ENABLE_OPENMP)

  IF (NOT OPENMP_FOUND)
    MESSAGE(FATAL_ERROR "No known compiler option for enabling OpenMP")
  ENDIF(NOT OPENMP_FOUND)

ENDIF(ENABLE_OPENMP)

SET(BUILD_PREFIX ${CMAKE_BINARY_DIR}/ep_build)

