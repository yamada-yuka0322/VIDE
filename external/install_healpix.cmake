execute_process(COMMAND ${MAKE} -C ${HEALPIX_DIR} hdrcopy)

set(bdir ${HEALPIX_DIR}/sampler)

file(MAKE_DIRECTORY ${DEST_DIR})
file(COPY ${bdir}/bin ${bdir}/lib ${bdir}/include DESTINATION ${DEST_DIR} NO_SOURCE_PERMISSIONS)
