SET(STAMP ${BUILD_PREFIX}/patch_applied)

SET(RESULT "not ok")

IF (EXISTS ${STAMP})
  FILE(READ ${STAMP} RESULT)
ENDIF (EXISTS ${STAMP})

IF(NOT "${RESULT}" STREQUAL "ok")
  EXECUTE_PROCESS(COMMAND patch -p0 -i ${PATCH_FILE}
                  WORKING_DIRECTORY ${SOURCE_PREFIX}
                  RESULT_VARIABLE okcode)
  IF(NOT "${okcode}" STREQUAL "0")
    FILE(WRITE ${STAMP} "not-applied")
    MESSAGE(FATAL_ERROR "Patch not applied")
  ELSE(NOT "${okcode}" STREQUAL "0")
    FILE(WRITE ${STAMP} "ok")
  ENDIF(NOT "${okcode}" STREQUAL "0")
ENDIF(NOT "${RESULT}" STREQUAL "ok")
