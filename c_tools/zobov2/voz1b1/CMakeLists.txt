add_executable(voz1b1_2 voz1b1.cpp vozutil.c voz_io.cpp)
target_link_libraries(voz1b1_2 ${QHULL_LIBRARY} ${COSMOTOOL_LIBRARY} ${MATH_LIB})

include_directories(${COSMOTOOL_INCLUDE_PATH} ${QHULL_INCLUDE_PATH})
