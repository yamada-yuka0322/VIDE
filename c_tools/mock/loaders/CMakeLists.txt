SET(SIMU_LOADERS_SRCS 
   simulation_loader.cpp ramses_loader.cpp 
   gadget_loader.cpp flash_loader.cpp)

IF (SDF_SUPPORT)
 add_definitions(-DSDF_SUPPORT)
 SET(SIMU_LOADERS_SRCS ${SIMU_LOADERS_SRCS} sdf_loader.cpp multidark_loader.cpp)
ENDIF (SDF_SUPPORT)


add_library(simu_loaders ${SIMU_LOADERS_SRCS})
