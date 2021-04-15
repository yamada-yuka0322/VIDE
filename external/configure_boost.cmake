MESSAGE(STATUS "Reconfiguring using compiler flags=\"${COMPILER_EXTRA_FLAGS}\" and linker flags=\"${LINKER_EXTRA_FLAGS}\"")
FILE(WRITE ${SRC_DIR}/tools/build/src/user-config.jam
  "using ${TOOLSET} : cmake : \"${COMPILER}\" : <cxxflags>\"${COMPILER_EXTRA_FLAGS}\" <linkflags>\"${LINKER_EXTRA_FLAGS}\" ;\n"
)

#FILE(APPEND ${SRC_DIR}/boost/regex/user.hpp
#  "#define BOOST_REGEX_MATCH_EXTRA 1\n"
#)

execute_process(
  COMMAND "${SRC_DIR}/bootstrap.sh" --prefix=${INSTALL_PATH} toolset=${TOOLSET}-cmake
  RESULT_VARIABLE okcode
)

IF (NOT "${okcode}" STREQUAL "0")
  MESSAGE(FATAL_ERROR "Cannot execute configure command")
ENDIF()
