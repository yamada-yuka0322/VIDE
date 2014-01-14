#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#ifdef PARALLEL
#include <mpi.h>
#endif

#include "header.h"

/* test file routines.
 */
void test(int argc, char **argv);
void chi2_grid(int argc, char **argv);
void fit_scale_bias(int argc, char **argv);
void aspen_breakout(void);
void populate_simulation_clf(void);

int main(int argc, char **argv)
{
  double s1;
  int i;

#ifdef PARALLEL
  printf("STARTING>>>\n");
  fflush(stdout);
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &ThisTask);
  MPI_Comm_size(MPI_COMM_WORLD, &NTask);
  printf("TASK %d reporting for duty.\n",ThisTask);
  fflush(stdout);
#endif

  OUTPUT=0;

  for(i=1;i<=99;++i)
    HOD.free[i]=0;
  wp.esys=0;

  Work.chi2=0;
  Work.imodel=1;

  USE_ERRORS = 0;
  ITRANS=4;
  HUBBLE=0.7;
  BEST_FIT = 0;
  HOD.M_sat_break = 1.0e14;
  HOD.alpha1 = 1.0;

  if(argc==1)
    endrun("./HOD.x hod.bat_file > output");

  read_parameter_file(argv[1]);

  /* If there's no cross-correlation function,
   * set the second number density equal to the first
   */
  if(!XCORR)
    GALAXY_DENSITY2 = GALAXY_DENSITY;

  /* Initialize the non-linear power spectrum.
   */
  nonlinear_sigmac(8.0);
  sigmac_interp(1.0E13);
  sigmac_radius_interp(1.0);

  if(argc>2)
    IDUM_MCMC=atoi(argv[2]);
  //  if(MCMC)m2n_mcmc();

  /* Get the galaxy bias factor
   */
  s1=qromo(func_galaxy_bias,log(HOD.M_low),log(HOD.M_max),midpnt);
  GALAXY_BIAS=s1/GALAXY_DENSITY;
  if(OUTPUT)
    fprintf(stdout,"Galaxy Bias bg= %f\n",GALAXY_BIAS);

  /* Get the galaxy satellite fraction
   */
  s1=qromo(func_satellite_density,log(HOD.M_low),log(HOD.M_max),midpnt)/
    GALAXY_DENSITY;
  if(OUTPUT)
    fprintf(stdout,"fsat %e\n",s1);

  /* Mean halo mass.
   */
  if(OUTPUT)
    fprintf(stdout,"M_eff %e\n",number_weighted_halo_mass());

  /* Set up BETA for wp integration.
   */
  BETA = pow(OMEGA_M,0.6)/GALAXY_BIAS;
  if(OUTPUT)
    printf("BETA = %f\n",BETA);

  /* Check for extra commands:
   * arg==999 goes to the test program, superceding tasks.
   * arg<0 supercedes the MCMC random number in the batfile.
   */
  if(argc>2)
    {
      if(atoi(argv[2])==999)
	test(argc,argv);
    }

  tasks(argc,argv);
}

