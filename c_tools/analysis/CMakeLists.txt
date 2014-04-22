include_directories(${CMAKE_CURRENT_BINARY_DIR})

SET(voidOverlap_SRCS voidOverlap.cpp)
add_genopt(voidOverlap_SRCS voidOverlap.ggo voidOverlap_conf STRUCTNAME voidOverlap_info)
add_executable(voidOverlap ${voidOverlap_SRCS})
target_link_libraries(voidOverlap  ${ZOB_LIBS})
