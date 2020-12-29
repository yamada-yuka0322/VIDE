SET(ENV{CFITSIO_LIBS} "-L${HEALPIX_INSTALL}/lib -lcfitsio")
SET(ENV{CFITSIO_CFLAGS} -I${HEALPIX_INSTALL}/include)

execute_process(
  COMMAND "${HEALPIX_DIR}/configure" --prefix=${HEALPIX_INSTALL} --libdir=${HEALPIX_INSTALL}/lib --disable-shared --with-pic --disable-maintainer-mode CC=${HEALPIX_CC} CXX=${HEALPIX_CXX}
    
    RESULT_VARIABLE okcode
)

IF (NOT "${okcode}" STREQUAL "0")
  MESSAGE(FATAL_ERROR "Configure failed")
ENDIF (NOT "${okcode}" STREQUAL "0")
