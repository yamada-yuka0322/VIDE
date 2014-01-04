#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "header.h"


/* This is the bias factor for halos. Mass (in h^{-1} M_sol) is the input.
 * This is not the scale-dependent bias factor, which can be calculated later.
 * 
 * There are many variants of the bias in this function. 
 * The one it is set up to use is from the Tinker etal (in prep, currently) 
 * bias paper calibrated from the SO halo catalogs from the Tinker et al mass functino
 * paper. The parameters of the bias functions are calibrated such that they will
 * work with any chosen halo overdensity.
 *
 * Other variants:
 *
 * The parameterization of b(M) is taken from Sheth Mo & Tormen (2001 MNRAS 232 1)
 * but the actual parameters have been replaced by the values of Appendix A in 
 * Tinker, Weinberg, Zheng, & Zehavi. 2005 Apj (M/L paper)
 *
 * See also Seljak & Warren (2004).
 * See also Sheth & Tormen (1999)
 * See also Mo & White (1996)
 */

double bias_from_file(double m, double r);

double bias(double mass)
{
  double rm,sig,k,neff,b,logk,khi,klo,phi,plo,nu,psp,x;
  static int flag=0;
  static double pnorm, prev_delta=-1;

  // variables for the SO(DELTA) bias functions
  static double bias_A, bias_a, bias_B, bias_b;

  /* Original SMT parameters */
  double a=0.707,bb=0.5,c=0.6;

  pnorm=SIGMA_8/sigmac(8.0);
  rm=pow(3.0*mass/(4.0*PI*OMEGA_M*RHO_CRIT),1.0/3.0);
  sig=pnorm*sigmac(rm);

  /* Tinker et al parameters */
  a=0.707;
  bb=0.35;
  c=0.8;

  /* Fitting to Mike Warren's simulations. */
  bb=0.28; 

  /* Use the globel parameters. */
  a=BIAS_A; bb=BIAS_B; c=BIAS_C;

  /* First normalize the power spectrum
   */
  pnorm=SIGMA_8/sigmac(8.0);
  rm=pow(3.0*mass/(4.0*PI*OMEGA_M*RHO_CRIT),1.0/3.0);
  sig=pnorm*sigmac(rm);



  /* This is from Tinker etal in prep for SO halos
   */
  if(DELTA_HALO != prev_delta)
    {
      bias_A = pow(fabs(log10(DELTA_HALO)-2.69),2)*0.16 + 0.785;
      bias_a = (log10(DELTA_HALO) - 2.28)*0.45;
      bias_B = 0.4*(log10(DELTA_HALO)-2.3)+0.7*log10(DELTA_HALO)*exp(-pow(4/log10(DELTA_HALO),4.0))+1.23;
      bias_b = 2.4;
      fprintf(stderr,"BIAS PARAMS: %f %f %f %f\n",bias_A, bias_a, bias_B, bias_b);
      prev_delta = DELTA_HALO;
    }
  
  a = pow(sig,-bias_a);
  b = 1 - bias_A*a/(a+1) + bias_B*pow(sig,-bias_b);

  return(b);

  /* This is the Seljak & Warren (2004) bias.
   * There's an alpha_s in the correction term, which here I set to zero. (RUNNING.)
   * See App. A of Tinker et al (M/L) for a comparison of this bias with above bias.
   */
  x=mass/MSTAR;
  b=(0.53+0.39*pow(x,0.45)+0.13/(40*x+1) + 5.0e-4*pow(x,1.5));
  b=b+log10(x)*(0.4*(OMEGA_M-0.3+SPECTRAL_INDX-1)+0.3*(SIGMA_8-0.9+HUBBLE-0.7)+0.8*0.0);
  return(b);

  /* This is the Sheth Mo Tormen bias.
   * (possible that parameter values have been changed to better fit simulations,
   * ie from Tinker etal 2005 ML paper).
   */
  nu=DELTA_CRIT/sig;
  b=1+1.0/(sqrt(a)*DELTA_CRIT)*(sqrt(a)*a*nu*nu + sqrt(a)*bb*pow(a*nu*nu,1-c) - 
				(pow(a*nu*nu,c)/(pow(a*nu*nu,c)+bb*(1-c)*(1-c/2.))));
  return(b);



  /* This is Sheth & Tormen
   */
  nu = DELTA_CRIT/sig;
  nu = nu*nu;
  return(1 + (0.707*nu - 1)/DELTA_CRIT + 2*0.3/DELTA_CRIT/(1+pow(0.707*nu,0.3)));
 
  /* This is the old Mo & White (1996) formula
   */
  return(1+DELTA_CRIT/sig/sig-1/DELTA_CRIT);
 

}

/* Just like the halo mass function, we'll set this up such that
 * you just need to interpolate.
 *
 * Now this has been changed to calculate the spatial-scale dependence
 * of the bias factor. If you don't want the scale dep. b, then just
 * input r<0.
 *
 * If the global flag LINEAR_PSP==0, uses the scale dependence calculated for
 * for halo bias relative to the non-linear matter \xi_m(r):
 *
 * f^2(r) = (1.0+xi*1.17)^1.49/(1.0+xi*0.69)^2.09  --> b(r) = b0*f(r)
 *
 * For LINEAR_PSP==1, use scale dependence determined for the linear P(k):
 *
 * f(r) = 1 + exp[-(r/A)^0.7] --> where A is a parameter that we're gonna have to 
 *                                determine in more detail, but now A \approx 1
 */
double bias_interp(double m, double r)
{
  static int flag=0,prev_cosmo=0;
  static double *x,*y,*y2, pnorm;
  int i,n=100;
  double dm,max=16.3,min=9,a,b,m1,m2,dm1,xi,power,rm,sig,b1,b2,mass,rvir,a1;

  if(!flag || RESET_COSMOLOGY!=prev_cosmo)
    {
      if(!ThisTask && OUTPUT)
	fprintf(stdout,"RESET: resetting bias for %f %f\n",OMEGA_M,SIGMA_8);
      if(!flag)
	{
	  x=dvector(1,n);
	  y=dvector(1,n);
	  y2=dvector(1,n);
	}
      flag=1;
      dm=(double)(max-min)/n;
      for(i=1;i<=n;++i)
	{
	  x[i]=pow(10.0,min+i*dm);
	  y[i]=log(bias(x[i]));
	  x[i] = log(x[i]);
	  //printf("BIAS %f %f\n",x[i]/2.302, y[i]/2.302);
	  continue;

	  // no longer need to do this part, since we're taking into account
	  // halo overdensity in the bias formula.
	  if(DELTA_HALO!=200)
	    {
	      x[i]=log(halo_mass_conversion2(x[i],halo_c200(x[i]),200.0,DELTA_HALO));
	    }
	  else
	    {
	      x[i]=log(x[i]);
	    }
	}
      spline(x,y,n,2.0E+30,2.0E+30,y2);
      prev_cosmo=RESET_COSMOLOGY;
      pnorm=SIGMA_8/sigmac(8.0);

    }


  m=log(m);
  splint(x,y,y2,n,m,&a);
  a = exp(a);

  // if we're using systematic errors in an MCMC, adjust parameter a1 (amplitude)
  if(USE_ERRORS)
    a *= M2N.bias_amp;

  // if no scale-dependence required, return the large-scale value
  if(r<0) return a;


  /* 
   * SCALE DEPENDENT BIAS!!!
   *----------------------------------------------------------------------
   */

  /* first, calculate the standard Zheng-like, mass-independent scale bias
   * based on the non-linear correlation function
   * (re-calibrated using the SO catalogs)
   */
  xi = xi_interp(r);  
  b = pow(1.0+xi*0.92,2.08)*pow(1.0+xi*0.74,-2.37);
  if(b<0.6)b=0.6;

  /* Now the mass-dependent term.
   * Where the mass-dependence comes in the form of the large-scale bias itself
   */
  /*
  sig = sigmac_radius_interp(r);  
  if(a<0)
    b *= pow(1 + 0.028*pow(a,1.53)*pow(sig,1.57),2.59)/
      pow(1 + 0.253*pow(a,1.24)*pow(sig,1.71),0.397);
  */

  // if we're using systematic errors in an MCMC, adjust parameter a1 (amplitude)
  if(USE_ERRORS) {
    b = (b-1)*M2N.scalebias_amp + 1;
    if(b<0)b=0;
  }
  return b*a;

}


/* This is the integrand which qromo or qtrap would call
 * to calculate the large-scale galaxy bias.
 * The integral is a number-weighted average of the halo
 * bias function, integrated over the halo mass function.
 */
double func_galaxy_bias(double m)
{
  m=exp(m);
  return(dndM_interp(m)*N_avg(m)*bias_interp(m,-1.)*m);
}

/*
 * This is to get the bias from a series of files.
 *
 */
double bias_from_file(double m, double r)
{
  FILE *fp;
  static int nfiles=10, *n, flag =1;
  static double *mass, **rad, **bias;
  float x1,x2,x3;
  char fname[1000];
  int i,j,i1,i2;
  double b;

  if(flag)
    {
      n = ivector(0,nfiles-1);
      mass = dvector(0,nfiles-1);
      rad = dmatrix(0,nfiles-1,1,30);
      bias = dmatrix(0,nfiles-1,1,30);

      fp = openfile("/home/tinker/cosmo/SDSS/zspace_analysis/SO_calibration/halofiles/sigma.dat");

      for(i=0;i<nfiles;++i)
	fscanf(fp,"%lf %f",&mass[i],&x1);
      for(i=0;i<nfiles;++i)
	mass[i] = log(mass[i]);
      fclose(fp);

      for(i=0;i<nfiles;++i)
	{
	  sprintf(fname,"/home/tinker/cosmo/SDSS/zspace_analysis/SO_calibration/halofiles/bias_nl.200.m%d",i);
	  fp = openfile(fname);
	  n[i] = filesize(fp);
	  for(j=1;j<=n[i];++j)
	    {
	      fscanf(fp,"%lf %lf %f %f %f",&rad[i][j],&bias[i][j],&x1,&x2,&x3);
	      rad[i][j] = log(rad[i][j]);
	    }
	  fclose(fp);
	}
      flag = 0;
    }

  r = log(r);
  if(m>mass[nfiles-1])return -1;

  for(i=1;i<nfiles;++i)
    if(mass[i]>=m)break;
  i2 = i;
  i1 = i-1;
  

  for(j=2;j<=n[i1];++j)
    if(rad[i1][j]>r)break;
  x1 = (bias[i1][j] - bias[i1][j-1])/(rad[i1][j] - rad[i1][j-1])*(r-rad[i1][j-1]) + bias[i1][j-1];

  for(j=2;j<=n[i2];++j)
    if(rad[i2][j]>r)break;
  x2 = (bias[i2][j] - bias[i2][j-1])/(rad[i2][j] - rad[i2][j-1])*(r-rad[i2][j-1]) + bias[i2][j-1];

  b = (x2 - x1)/(mass[i2] - mass[i1])*(m - mass[i1]) + x1;
  
  //  if(r>log(4))
  //  printf("%e %e %f %f %f %f %f\n",m,r,b,x2,x1,mass[i2],mass[i1]);
  return b;

}
