#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "header.h"

double collapse_redshift(double z);
double cvir_pnorm_g1,
  cvir_sigma_g1;

/* This calculates and tabulates the halo concentrations
 * as a function of halo mass. Uses the "Bullock model", 
 * described in a little more detail below.
 */

double halo_concentration(double m)
{
  static int flag=1,n,prev_cosmo=0;
  static double *x,*y,*y2;
  int i;
  float x1,x2,cfac;
  double a,dm,x3,x4;
  FILE *fp;
  char fname[1000];

  if(flag || RESET_COSMOLOGY!=prev_cosmo)
    {
      MSTAR = mstar();
      if(OUTPUT)
	fprintf(stdout,"Calc cvir with DELTA_HALO= %f\n",DELTA_HALO);
      prev_cosmo=RESET_COSMOLOGY;
      n=50;
      if(flag)
	{
	  x=dvector(1,n);
	  y=dvector(1,n);
	  y2=dvector(1,n);
	}
      cvir_pnorm_g1=SIGMA_8/sigmac(8.0);

      flag=0;
      
      dm=(log(HOD.M_max)-log(1.0E8))/(n-1);
      for(i=1;i<=n;++i)
	{
	  x[i]=exp((i-1)*dm+log(1.0E8));
	  y[i]=cvir_model(x[i]);
	  x[i]=log(halo_mass_conversion(x[i],&y[i],DELTA_HALO));
	  y[i]=log(y[i]);
	  //y[i] = log(10);
	}
      spline(x,y,n,1.0E+30,1.0E+30,y2);
    }
  cfac = 0.13*log10(m/1.0e12) + 0.125;
  cfac = 1;
  if(m>80*MSTAR) m=80*MSTAR;
  m=log(m);
  splint(x,y,y2,n,m,&a);
  return(exp(a)*cfac);
  
}


/* This is the cvir model, which should reproduce the results
 * of cvir3.f, which can be obtained off James Bullock's web site.
 *
 * Basic idea; For a halo of mass M, find the collapse redshift
 * of that halo be finding the redshift at which rms mass variance
 * sigma(Mf) = DELTA_CRIT (using linear growth factor), where Mf is some
 * fraction of the redshift zero halo mass.
 *
 * cvir = k*(1+z_coll(M*f))
 *
 * Model params:
 *  - k = 3.12    - fitting parameter (taken from cvir3.f)
 *  - f = 0.001   - fraction of halo mass that collapsed at z_coll
 */

double cvir_model(double mass)
{
  double cvir,z_coll,zbrent(),rad;
  double k=3.12;
  double f=0.001;
 
  rad=pow(3.0*f*mass/(4.0*PI*OMEGA_M*RHO_CRIT),1.0/3.0);
  cvir_sigma_g1=cvir_pnorm_g1*sigmac(rad);

  if(collapse_redshift(0.0)<0)return(k);
  z_coll=zbrent(collapse_redshift,0.0,200.0,1.0E-3);
  cvir=k*(1+z_coll);
  return(cvir);
}

/* This is the input function to find where sig(M)*D(z)-DELTA_CRIT = 0
 * which is the redshift at which a halo of mass M collapsed.
 */
double collapse_redshift(double z)
{
  double D;
  D=growthfactor(z);
  return(cvir_sigma_g1*D-DELTA_CRIT);
}

/* Some quantities are specific to an overdensity of 200 (i.e., the Jenkins mass
 * function and the halo bias relation of Tinker et al. )
 * Specfically, these are actually for FOF 0.2 halos, for which 200 is the
 * current best approximation. (Although Warren et al quotes 250, which is the most recent
 * value.)
 *
 * Therefore, the halo concentrations for Delta=200 halos need to be easily 
 * accessible for halo mass conversion of these quantities. The values are 
 * tabulates here, as opposed to the above routine which tabulates the 
 * concentrations for a user-specified overdensity.
 */

double halo_c200(double m)
{
  static int flag=1,n,prev_cosmo=0;
  static double *x,*y,*y2;
  int i;
  float x1,x2;
  double a,dm,x3,x4;
  FILE *fp;
  char fname[1000];

  if(flag || RESET_COSMOLOGY!=prev_cosmo)
    {
      if(OUTPUT)
	fprintf(stdout,"Calc cvir with DELTA_HALO= %f\n",200.0);
      prev_cosmo=RESET_COSMOLOGY;
      n=50;
      if(flag)
	{
	  x=dvector(1,n);
	  y=dvector(1,n);
	  y2=dvector(1,n);
	}
      cvir_pnorm_g1=SIGMA_8/sigmac(8.0);

      flag=0;
      
      dm=(log(HOD.M_max)-log(1.0E8))/(n-1);
      for(i=1;i<=n;++i)
	{
	  x[i]=exp((i-1)*dm+log(1.0E8));
	  y[i]=cvir_model(x[i]);
	  x[i]=log(halo_mass_conversion(x[i],&y[i],200.0));
	  y[i]=log(y[i]);
	}
    }
  m=log(m);
  splint(x,y,y2,n,m,&a);
  return(exp(a));
  
}
