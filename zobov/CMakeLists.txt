
add_executable(voz1b1 voz1b1.c readfiles.c vozutil.c voz.h)
target_link_libraries(voz1b1 ${QHULL_LIBRARY} ${MATH_LIB})

add_executable(2Dvoz1b1 2Dvoz1b1.c readfiles.c vozutil.c voz.h)
target_link_libraries(2Dvoz1b1 ${QHULL_LIBRARY} ${MATH_LIB})


#add_executable(jozov jozov.c findrtop.c)
#target_link_libraries(jozov ${MATH_LIB})

add_executable(vozinit vozinit.c readfiles.c)
target_link_libraries(vozinit ${MATH_LIB})

add_executable(2Dvozinit 2Dvozinit.c readfiles.c)
target_link_libraries(2Dvozinit ${MATH_LIB})


add_executable(voztie voztie.c readfiles.c)
target_link_libraries(voztie ${MATH_LIB})

add_executable(2Dvoztie 2Dvoztie.c readfiles.c)
target_link_libraries(2Dvoztie ${MATH_LIB})

install(TARGETS voz1b1 vozinit voztie 2Dvoz1b1 2Dvozinit 2Dvoztie DESTINATION ${VIDE_BIN})
