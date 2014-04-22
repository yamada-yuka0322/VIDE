include_directories(${CMAKE_CURRENT_BINARY_DIR})

SET(stackVoids_SRCS stackVoids.cpp)
add_genopt(stackVoids_SRCS stackVoids.ggo stackVoids_conf STRUCTNAME stackVoids_info)
add_executable(stackVoids ${stackVoids_SRCS})
target_link_libraries(stackVoids  ${ZOB_LIBS})

SET(stackVoidsZero_SRCS stackVoidsZero.cpp)
add_genopt(stackVoidsZero_SRCS stackVoidsZero.ggo stackVoidsZero_conf STRUCTNAME stackVoidsZero_info)
add_executable(stackVoidsZero ${stackVoidsZero_SRCS})
target_link_libraries(stackVoidsZero ${ZOB_LIBS})

set(computeVelocityProfile_SRCS computeVelocityProfile.cpp)
add_genopt(computeVelocityProfile_SRCS computeVelocityProfile.ggo computeVelocityProfile_conf STRUCT Params)
add_executable(computeVelocityProfile ${computeVelocityProfile_SRCS})
target_link_libraries(computeVelocityProfile  ${ZOB_LIBS})


set(stackDensityField_SRCS stackDensityField.cpp)
add_genopt(stackDensityField_SRCS stackDensityField.ggo stackDensityField_conf STRUCT PARAMS)
add_executable(stackDensityField ${stackDensityField_SRCS})
target_link_libraries(stackDensityField ${ZOB_LIBS})



set(stackVelocityField_SRCS stackVelocityField.cpp)
add_genopt(stackVelocityField_SRCS stackVelocityField.ggo stackVelocityField_conf STRUCT PARAMS)
add_executable(stackVelocityField ${stackVelocityField_SRCS})
target_link_libraries(stackVelocityField ${ZOB_LIBS})

SET(pruneVoids_SRCS pruneVoids.cpp)
add_genopt(pruneVoids_SRCS pruneVoids.ggo pruneVoids_conf STRUCTNAME pruneVoids_info)
add_executable(pruneVoids ${pruneVoids_SRCS})
target_link_libraries(pruneVoids ${ZOB_LIBS})

SET(makeAHFOutput_SRCS makeAHFOutput.cpp)
add_genopt(makeAHFOutput_SRCS makeAHFOutput.ggo makeAHFOutput_conf STRUCTNAME makeAHFOutput_info)
add_executable(makeAHFOutput ${makeAHFOutput_SRCS})
target_link_libraries(makeAHFOutput ${ZOB_LIBS})
