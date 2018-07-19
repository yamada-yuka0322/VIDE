execute_process(COMMAND "${PYTHON_EXECUTABLE}" "-c"
	"import distutils.sysconfig as cs; import os; import sys; v=cs.get_config_vars(); print(os.path.join(v['LIBDIR'],v['LDLIBRARY'])); sys.exit(0)"
    RESULT_VARIABLE _PYLIB_SEARCH_SUCCESS
    OUTPUT_VARIABLE _PYLIB_VALUES_OUTPUT
    ERROR_VARIABLE _PYLIB_ERROR_VALUE
    OUTPUT_STRIP_TRAILING_WHITESPACE)

message(${_PYLIB_SEARCH_SUCCESS})

execute_process(COMMAND "${PYTHON_EXECUTABLE}" "-c"
    "import distutils.sysconfig as cs;  import os; v=cs.get_config_vars(); print(v['INCLUDEPY']);"
    RESULT_VARIABLE _PYINC_SEARCH_SUCCESS
    OUTPUT_VARIABLE _PYINC_VALUES_OUTPUT
    ERROR_VARIABLE _PYINC_ERROR_VALUE
    OUTPUT_STRIP_TRAILING_WHITESPACE)


if(NOT _PYLIB_SEARCH_SUCCESS MATCHES 0)
    message(FATAL_ERROR
           "PyLib search failure:\n${_PYLIB_ERROR_VALUE}")
    return()
endif()

if(NOT _PYINC_SEARCH_SUCCESS MATCHES 0)
    message(FATAL_ERROR
           "PyInc search failure:\n${_PYINC_ERROR_VALUE}")
    return()
endif()


set(PYTHON_LIBRARY ${_PYLIB_VALUES_OUTPUT} CACHE PATH "Python runtime library path")
set(PYTHON_INCLUDE_PATH ${_PYINC_VALUES_OUTPUT} CACHE PATH "Python runtime include path")
