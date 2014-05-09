cmake_minimum_required(VERSION 2.8)

project(zobovPerso)

include(CheckCCompilerFlag)
include(ExternalProject)

check_c_compiler_flag(-march=native SUPPORT_ARCH_NATIVE )



find_library(MATH_LIB m)

macro(add_genopt _sourcelist _ggofile _basefile)

  unset(_structname)
  unset(_funcname)

  if(NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/${_ggofile})
    set(_ggofile2 ${CMAKE_CURRENT_SOURCE_DIR}/${_ggofile})
  else(NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/${_ggofile})
    set(_ggofile2 ${CMAKE_CURRENT_BINARY_DIR}/${_ggofile})
  endif(NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/${_ggofile})

  foreach(arg ${ARGN})
    if ("x${arg}" MATCHES "^x(STRUCTNAME|FUNCNAME)$")
      SET(doing "${arg}")
    elseif(doing STREQUAL "STRUCTNAME")
      SET(_structname ${arg})
    elseif(doing STREQUAL "FUNCNAME")
      SET(_funcname ${arg})
    endif()
  endforeach(arg ${ARGN})

  if(NOT DEFINED _structname)
    set(_structname ${_basefile})
  endif(NOT DEFINED _structname)

  if(NOT DEFINED _funcname)
    set(_funcname ${_basefile})
  endif(NOT DEFINED _funcname)
  
  set(_cfile ${CMAKE_CURRENT_BINARY_DIR}/${_basefile}.c)
  set(_hfile ${CMAKE_CURRENT_BINARY_DIR}/${_basefile}.h)

  add_custom_command(
    OUTPUT ${_cfile} ${_hfile}
    COMMAND ${GENGETOPT} -i ${_ggofile2} -f ${_funcname} -a ${_structname} zobovConf_info -F ${_basefile} -C
    DEPENDS ${_ggofile2}
    )

  set(${_sourcelist} ${_cfile} ${${_sourcelist}})  

endmacro(add_genopt)

macro(configure_exec _source _destdir _destfile)

  SET(TMP_CONFIGURE_DIR ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY})
  configure_file(${_source} ${TMP_CONFIGURE_DIR}/${_destfile} @ONLY)
  file(COPY ${TMP_CONFIGURE_DIR}/${_destfile} DESTINATION ${_destdir} 
       FILE_PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
 
endmacro(configure_exec)

include(${CMAKE_SOURCE_DIR}/external/external_build.cmake)
include(${CMAKE_SOURCE_DIR}/external/external_python_build.cmake)

SET(PYTHONPATH ${PYTHON_LOCAL_SITE_PACKAGE})
configure_exec(${CMAKE_SOURCE_DIR}/run_python.sh.in ${CMAKE_BINARY_DIR} run_python.sh)
configure_exec(${CMAKE_SOURCE_DIR}/run_env.sh.in  ${CMAKE_BINARY_DIR} run_env.sh)
configure_exec(${CMAKE_SOURCE_DIR}/python_tools/pipeline_source/prepareInputs.in.py
${CMAKE_BINARY_DIR}/pipeline prepareInputs.py)
#configure_exec(${CMAKE_SOURCE_DIR}/python_tools/pipeline_source/applyMaskToMock.in.py
${CMAKE_BINARY_DIR}/pipeline applyMaskToMock.py)
#configure_exec(${CMAKE_SOURCE_DIR}/python_tools/pipeline_source/buildSkyProjections.in.py
${CMAKE_BINARY_DIR}/pipeline/miscTools buildSkyProjections.py)

SET(python_build_environment
  ${CMAKE_COMMAND} -DPYTHON_LOCAL_SITE_PACKAGE=${PYTHON_LOCAL_SITE_PACKAGE} -DVOID_GSL=${CMAKE_BINARY_DIR}/ext_build/gsl -DPYTHON_EXECUTABLE=${PYTHON_EXECUTABLE} -DTARGET_PATH=${CMAKE_BINARY_DIR}/ext_build/python -P)

add_custom_target(python_pipeline ALL
  COMMAND ${python_build_environment} ${CMAKE_SOURCE_DIR}/external/python_build.cmake
  COMMAND ${python_build_environment} ${CMAKE_SOURCE_DIR}/external/python_install.cmake
  DEPENDS gsl cython netcdf4-python ${PYTHON_AUXILIARY_DEPEND}
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/python_tools
)


subdirs(zobov c_tools)
