include_directories(${CMAKE_CURRENT_BINARY_DIR})

SET(hod_SRCS header.c main.c utility.c sigmac.c transfnc.c transfunc_file.c 
  nonlinear_power_spectrum.c least_squares.c hod_functions.c 
   xi_matter.c one_halo_rspace.c two_halo_rspace.c 
  input_params.c dFdx.c mstar.c halo_concentration.c growthfactor.c 
  halo_mass_conversion.c nfw_transform.c 
  pair_density.c tasks.c wp_minimization.c 
  kaiser_distortions.c 
  tf_eisenstein_hu.c  jeans.c
  populate_simulation.c test.c
  dark_matter_statistics.c 
  meshlink2.c nbrsfind2.c i3tensor_2.c 
  mcmc.c halo_mass_function.c halo_bias.c
  nrutil.c qromo.c midpnt.c midinf.c polint.c splint.c spline.c 
  zbrent.c qtrap.c trapzd.c cisi.c complex.c amoeba.c amotry.c 
  gaussj.c powell.c linmin.c f1dim.c mnbrak.c brent.c gasdev.c 
  ran1.c jacobi.c splin2.c splie2.c ran2.c sort2.c
)
add_executable(hod ${hod_SRCS})
target_link_libraries(hod  -lm)
