#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#ifdef PARALLEL
#include <mpi.h>
#endif
#include "header.h"

/***********************************************************************
 * This is the main function for controlling the activities of the code.
 * If you ask for any sort of minimization or MCMC analysis, the code is
 * redirected to another routine, but once finished will complete the rest
 * of the tasks checked off (so if you want to minimize and then print out
 * related statistics, the code will do that automatically if you ask it).
 *
 *
 * Each task generates an output file using the "root_filename" in the hod.bat
 * file and tacking on an extension at the end:
 *
 *  [filename].r_space  --> real space correlation function
 *  [filename].HOD      --> prints out the mean N(M)
 *
 * See the full list below.
 */


void tasks(int argc, char **argv)
{
  int i,j,nsize,ix,iy,nr,i1,i2,n;
  float x1,x2,x3,x4,x5,x6,err,npairs;
  double r,rs,rp,dx1,dx2,dx3,dx4,dx5,**xi2d,*xx2d,*yy2d,**xi2d_data,**xi2d_kaiser,
    xi0_m,xi2_m,xi0_k,xi2_k,xi0_d,xi2_d,xi_r,rlo,rhi,rr[50],rhalf[50],dr,
    rphi,rshi,rslo,rplo,dlogm,m,sig;
  FILE *fp,*fp2,*fp1;
  char fname[100];

  


  /* This is for chi^2 minimization of data for the projected correlation
   * function.
   */
  if(Task.wp_minimize && !HOD.color)
    wp_minimization(argv[1]);

  /* This is for Monte-Carlo Markov Chain exploration of the posterior
   * distribution of the parameters, either real-space or redshift-space,
   * depending on what MCMC is set to.
   */
  if(Task.MCMC)
    mcmc_minimization();

  /* This is to output the shape of the mean occupation function and the 
   * scatter about the mean. 
   * File columns are:
   *  1 - halo mass (M_sol/h)
   *  2 - N_cen(M)
   *  3 - N_sat(M)
   *  4 - N_tot(M)
   *  5 - <N(N-1)> 
   *  6 - ratio of scatter to Poisson scatter (often called "alpha")
   */
  if(Task.HOD || Task.All)
    {
      sprintf(fname,"%s.HOD",Task.root_filename);
      fp=fopen(fname,"w");
      dlogm=(log(HOD.M_max)-log(HOD.M_low))/99;
      for(i=1;i<=100;++i)
	{
	  m=exp((i-1)*dlogm)*HOD.M_low;
	  sig = N_sat(m)*N_sat(m) + N_cen(m)*(1-N_cen(m));
	  fprintf(fp,"%e %e %e %e %e %e\n",
		  m,N_cen(m),N_sat(m),N_avg(m),sig,sig/(N_avg(m)*N_avg(m)));
	}
      fclose(fp);
    }

  /* This sets in motion the calculation and tabulation of the real-space
   * correlation function and outputs it to a file.
   * File columns are:
   *  1 - radius (Mpc/h)
   *  2 - one-halo term (real-space)
   *  3 - two-halo term (real-space)
   *  4 - full xi(r)
   *  5 - projected correlation function (1st column is now r_p)
   *  6 - projected correlation function without z-space correction
   */
  if(Task.real_space_xi || Task.All)
    {
      fprintf(stderr,"\n\nCALCULATING REAL-SPACE CORRELATION FUNCTION.\n");
      fprintf(stderr,    "--------------------------------------------\n\n");

      sprintf(fname,"%s.r_space",Task.root_filename);
      fp=fopen(fname,"w");
      dr=(log(70.0)-log(0.01))/49.0;
      for(i=0;i<50;++i)
	{
	  r=exp(i*dr+log(0.01));
	  x1=one_halo_real_space(r);
	  x2=two_halo_real_space(r);
	  x3=projected_xi(r);
	  x4 = projected_xi_rspace(r);
	  fprintf(fp,"%f %e %e %e %e %e\n",r,x1,x2,x1+x2,x3,x4);
	  fflush(fp);
	}
      fclose(fp);
    }

  /* This takes a halofile from a simulation and populates the halos
   * with galaxies according the HOD specified in the batch file.
   */ 
  if(Task.populate_sim==1)
    populate_simulation();

  if(!ThisTask)
    OUTPUT=1;

  /* Output the linear and non-linear dark matter power spectrum.
   * Non-linear P(k) calculated using Smith et al.
   *
   * Format of file [root].matter_pk
   *   1 - k [h/Mpc]
   *   2 - linear P(k) [Mpc/h]^3
   *   3 - non-linear P(k) [Mpc/h]^3
   */
  if(Task.matter_pk)
    output_matter_power_spectrum();

  /* Output the linear and non-linear dark matter power spectrum.
   * Non-linear xi(r) is Fourier transform of Smith et al (above)
   *
   * Format of file [root].matter_pk
   *   1 - r [Mpc/h]
   *   2 - linear xi(r)
   *   3 - non-linear xi(r)
   */
  if(Task.matter_xi)
    output_matter_correlation_function();

  /* Output the matter variance as a function of scale.
   *
   * Format of file [root].sigma_r 
   *  1 - r [Mpc/h]
   *  2 - sigma(r) [linear]
   *  3 - sigma(r) [non-linear, using Smith]
   *  4 - mass [M_sol/h] mass = (4/3)*PI*r^3*rho_bar
   */
  if(Task.sigma_r)
    output_matter_variance();

  /* Output the halo concentrations using Bullock et al (2001) model
   * 
   * Format of file [root].civr
   *  1 - mass [Msol/h] --> mass specified by DELTA_HALO (input file).
   *  2 - halo concentration. (for the specified halo definition)
   */
  if(Task.cvir)
    output_halo_concentrations();

  /* Output the halo mass function using the results of Tinker et al (2008)
   *
   * Format of the file [root].massfunc
   *  1 - mass [Msol/h]
   *  2 - dn/dM [Mph/c]^-3 (the differential mass function).
   */
  if(Task.dndM)
    output_halo_mass_function();

  endrun("finished with tasks");
}
