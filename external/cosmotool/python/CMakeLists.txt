set(CMAKE_SHARED_MODULE_PREFIX)


set(PYTHON_INCLUDES ${NUMPY_INCLUDE_DIRS} ${PYTHON_INCLUDE_PATH} ${CMAKE_SOURCE_DIR}/python)
include_directories(${CMAKE_SOURCE_DIR}/src ${CMAKE_BINARY_DIR}/src)

IF(CYTHON)
  add_custom_command(
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/_cosmotool.cpp
    COMMAND ${CYTHON} --cplus -o ${CMAKE_CURRENT_BINARY_DIR}/_cosmotool.cpp ${CMAKE_CURRENT_SOURCE_DIR}/_cosmotool.pyx
    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/_cosmotool.pyx) 

  add_custom_command(
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/_cosmo_power.cpp
    COMMAND ${CYTHON} --cplus -o ${CMAKE_CURRENT_BINARY_DIR}/_cosmo_power.cpp ${CMAKE_CURRENT_SOURCE_DIR}/_cosmo_power.pyx
    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/_cosmo_power.pyx) 

  add_custom_command(
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/_fast_interp.cpp
    COMMAND ${CYTHON} --cplus -o ${CMAKE_CURRENT_BINARY_DIR}/_fast_interp.cpp ${CMAKE_CURRENT_SOURCE_DIR}/_fast_interp.pyx
    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/_fast_interp.pyx) 

  add_custom_command(
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/_cosmo_cic.cpp
    COMMAND ${CYTHON} --cplus -o ${CMAKE_CURRENT_BINARY_DIR}/_cosmo_cic.cpp ${CMAKE_CURRENT_SOURCE_DIR}/_cosmo_cic.pyx
    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/_cosmo_cic.pyx) 

  add_custom_command(
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/_project.cpp
    COMMAND ${CYTHON} --cplus -o ${CMAKE_CURRENT_BINARY_DIR}/_project.cpp ${CMAKE_CURRENT_SOURCE_DIR}/_project.pyx
    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/_project.pyx ${CMAKE_CURRENT_SOURCE_DIR}/project_tool.hpp ) 

ENDIF(CYTHON)


add_library(_cosmotool MODULE ${CMAKE_CURRENT_BINARY_DIR}/_cosmotool.cpp)
add_library(_cosmo_power MODULE ${CMAKE_CURRENT_BINARY_DIR}/_cosmo_power.cpp)
add_library(_cosmo_cic MODULE ${CMAKE_CURRENT_BINARY_DIR}/_cosmo_cic.cpp)
add_library(_fast_interp MODULE ${CMAKE_CURRENT_BINARY_DIR}/_fast_interp.cpp)
add_library(_project MODULE ${CMAKE_CURRENT_BINARY_DIR}/_project.cpp)
target_include_directories(_cosmotool PRIVATE ${PYTHON_INCLUDES}) 
target_include_directories(_cosmo_power PRIVATE ${PYTHON_INCLUDES}) 
target_include_directories(_cosmo_cic PRIVATE ${PYTHON_INCLUDES}) 
target_include_directories(_fast_interp PRIVATE ${PYTHON_INCLUDES}) 
target_include_directories(_project PRIVATE ${PYTHON_INCLUDES}) 


SET(CMAKE_MODULE_LINKER_FLAGS "${CMAKE_MODULE_LINKER_FLAGS} -Bsymbolic-functions")
if(APPLE)
  set(CMAKE_MODULE_LINKER_FLAGS "-undefined dynamic_lookup")
endif()

target_link_libraries(_cosmotool ${CosmoTool_local} ${PYTHON_LIBRARIES} ${GSL_LIBRARIES})
target_link_libraries(_cosmo_power ${CosmoTool_local} ${PYTHON_LIBRARIES} ${GSL_LIBRARIES})
target_link_libraries(_cosmo_cic ${CosmoTool_local} ${PYTHON_LIBRARIES} ${GSL_LIBRARIES})
target_link_libraries(_project ${PYTHON_LIBRARIES})
target_link_libraries(_fast_interp ${CosmoTool_local} ${PYTHON_LIBRARIES})

SET(ct_TARGETS   _cosmotool _project _cosmo_power _cosmo_cic _fast_interp )

if (Boost_FOUND)
  message(STATUS "Building bispectrum support (path = ${Boost_INCLUDE_DIRS})")
  include_directories(${Boost_INCLUDE_DIRS})
  add_library(_cosmo_bispectrum MODULE _cosmo_bispectrum.cpp)
  target_link_libraries(_cosmo_bispectrum ${MATH_LIBRARY})
  if(ENABLE_OPENMP)
    set_target_properties(_cosmo_bispectrum PROPERTIES COMPILE_FLAGS "${OpenMP_CXX_FLAGS}" LINK_FLAGS "${OpenMP_CXX_FLAGS}")
  endif()
  if (Boost_DEP)
    add_dependencies(_cosmo_bispectrum ${Boost_DEP})
  endif()
  SET(ct_TARGETS ${ct_TARGETS} _cosmo_bispectrum)
endif()

# Discover where to put packages
if (NOT PYTHON_SITE_PACKAGES)
  execute_process (COMMAND ${PYTHON_EXECUTABLE} -c "from distutils.sysconfig import get_python_lib; print(get_python_lib())" OUTPUT_VARIABLE internal_PYTHON_SITE_PACKAGES OUTPUT_STRIP_TRAILING_WHITESPACE)
  SET(SYSTEM_PYTHON_SITE_PACKAGES ${internal_PYTHON_SITE_PACKAGES} CACHE PATH "Path to the target system-wide site-package where to install python modules")

  execute_process (COMMAND ${PYTHON_EXECUTABLE} -c "from site import USER_SITE; print(USER_SITE)" OUTPUT_VARIABLE internal_PYTHON_SITE_PACKAGES OUTPUT_STRIP_TRAILING_WHITESPACE)
  SET(USER_PYTHON_SITE_PACKAGES ${internal_PYTHON_SITE_PACKAGES} CACHE PATH "Path to the target user site-package where to install python modules")

  mark_as_advanced(USER_PYTHON_SITE_PACKAGES SYSTEM_PYTHON_SITE_PACKAGES)
endif (NOT PYTHON_SITE_PACKAGES)

message(STATUS "System python site: ${SYSTEM_PYTHON_SITE_PACKAGES}")
message(STATUS "User python site: ${USER_PYTHON_SITE_PACKAGES}")

OPTION(INSTALL_PYTHON_LOCAL OFF)

IF (NOT INSTALL_PYTHON_LOCAL)
  SET(PYTHON_SITE_PACKAGES ${SYSTEM_PYTHON_SITE_PACKAGES})
ELSE (NOT INSTALL_PYTHON_LOCAL)
  SET(PYTHON_SITE_PACKAGES ${USER_PYTHON_SITE_PACKAGES})
ENDIF(NOT INSTALL_PYTHON_LOCAL)
cmessage(STATUS "Python install location: ${PYTHON_SITE_PACKAGES}")


if (WIN32 AND NOT CYGWIN)
  SET_TARGET_PROPERTIES(_cosmotool PROPERTIES SUFFIX ".pyd")
endif (WIN32 AND NOT CYGWIN)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cosmotool/config.py.in ${CMAKE_CURRENT_BINARY_DIR}/cosmotool/config.py @ONLY)

INSTALL(TARGETS 
  ${ct_TARGETS}
  LIBRARY DESTINATION ${PYTHON_SITE_PACKAGES}/cosmotool
)

INSTALL(DIRECTORY cosmotool ${CMAKE_CURRENT_BINARY_DIR}/cosmotool DESTINATION ${PYTHON_SITE_PACKAGES} 
        FILES_MATCHING PATTERN "*.py")
