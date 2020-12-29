function(omp_add_flags TARGET LANG)
   if(OPENMP_FOUND)
     if(NOT LANG)
       get_target_property(LANG ${TARGET} LINKER_LANGUAGE)
     endif()
     if(LANG MATCHES "C(XX)?")
       set_property(TARGET ${TARGET} APPEND
                    PROPERTY COMPILE_FLAGS ${OpenMP_${LANG}_FLAGS})
       set_property(TARGET ${TARGET} APPEND
                    PROPERTY LINK_FLAGS ${OpenMP_${LANG}_FLAGS})
     else()
       message(WARNING "omp_add_flags: target '${TARGET}'"
                      " link language '${LANG}' is not supported")
    endif()
  endif()
endfunction()

add_executable(jozov2 jozov2.cpp jozov2_io.cpp jozov2_zones.cpp jozov2_watershed.cpp findrtop.c)

if (ENABLE_OPENMP)
  Find_Package(OpenMP)
  omp_add_flags(jozov2 CXX)
  add_definitions(-DOPENMP)
ENDIF(ENABLE_OPENMP)

install(TARGETS jozov2 DESTINATION ${VIDE_BIN})
