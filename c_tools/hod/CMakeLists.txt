include_directories(${CMAKE_CURRENT_BINARY_DIR})

SET(hod_SRCS header.c main.c utility.c sigmac.c
  hod_functions.c 
  halo_bias.c
  transfnc.c
  nonlinear_power_spectrum.c
  halo_mass_function.c
  tf_eisenstein_hu.c
  transfunc_file.c
  halo_concentration.c
  mstar.c
  xi_matter.c
  halo_mass_conversion.c
  least_squares.c
  jeans.c
  growthfactor.c
  input_params.c dFdx.c
  tasks.c 
  populate_simulation.c
  meshlink2.c nbrsfind2.c i3tensor_2.c 
  nrutil.c qromo.c midpnt.c midinf.c polint.c splint.c spline.c 
  zbrent.c qtrap.c trapzd.c cisi.c complex.c amoeba.c amotry.c 
  gaussj.c powell.c linmin.c f1dim.c mnbrak.c brent.c gasdev.c 
  ran1.c jacobi.c splin2.c splie2.c ran2.c sort2.c
)
add_executable(hod ${hod_SRCS})
target_link_libraries(hod  -lm)
install(TARGETS hod DESTINATION ${VIDE_BIN})
