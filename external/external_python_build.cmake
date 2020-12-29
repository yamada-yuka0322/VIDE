INCLUDE(FindPythonInterp)

# Discover where to put packages
if (NOT PYTHON_SITE_PACKAGES)
  execute_process (
     COMMAND ${PYTHON_EXECUTABLE} -c "from distutils.sysconfig import get_python_lib; print(get_python_lib())"
     OUTPUT_VARIABLE internal_PYTHON_SITE_PACKAGES
     OUTPUT_STRIP_TRAILING_WHITESPACE
  )
  SET(SYSTEM_PYTHON_SITE_PACKAGES
     ${internal_PYTHON_SITE_PACKAGES}
     CACHE PATH "Path to the target system-wide site-package where to install python modules")

  execute_process (
     COMMAND ${PYTHON_EXECUTABLE} -c "from site import USER_SITE; print(USER_SITE)"
     OUTPUT_VARIABLE internal_PYTHON_SITE_PACKAGES
     OUTPUT_STRIP_TRAILING_WHITESPACE
  )
  SET(USER_PYTHON_SITE_PACKAGES
     ${internal_PYTHON_SITE_PACKAGES}
     CACHE PATH "Path to the target user site-package where to install python modules")

  mark_as_advanced(USER_PYTHON_SITE_PACKAGES SYSTEM_PYTHON_SITE_PACKAGES)
  message(STATUS "System python site: ${SYSTEM_PYTHON_SITE_PACKAGES}")
  message(STATUS "User python site: ${USER_PYTHON_SITE_PACKAGES}")
endif (NOT PYTHON_SITE_PACKAGES)

