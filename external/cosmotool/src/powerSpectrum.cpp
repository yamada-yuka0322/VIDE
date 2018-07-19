/*+
This is CosmoTool (./src/powerSpectrum.cpp) -- Copyright (C) Guilhem Lavaux (2007-2014)

guilhem.lavaux@gmail.com

This software is a computer program whose purpose is to provide a toolbox for cosmological
data analysis (e.g. filters, generalized Fourier transforms, power spectra, ...)

This software is governed by the CeCILL license under French law and
abiding by the rules of distribution of free software.  You can  use, 
modify and/ or redistribute the software under the terms of the CeCILL
license as circulated by CEA, CNRS and INRIA at the following URL
"http://www.cecill.info". 

As a counterpart to the access to the source code and  rights to copy,
modify and redistribute granted by the license, users are provided only
with a limited warranty  and the software's author,  the holder of the
economic rights,  and the successive licensors  have only  limited
liability. 

In this respect, the user's attention is drawn to the risks associated
with loading,  using,  modifying and/or developing or reproducing the
software by the user in light of its specific status of free software,
that may mean  that it is complicated to manipulate,  and  that  also
therefore means  that it is reserved for developers  and  experienced
professionals having in-depth computer knowledge. Users are therefore
encouraged to load and test the software's suitability as regards their
requirements in conditions enabling the security of their systems and/or 
data to be ensured and,  more generally, to use and operate it in the 
same conditions as regards security. 

The fact that you are presently reading this means that you have had
knowledge of the CeCILL license and that you accept its terms.
+*/

#include <cassert>
#include <vector>
#include <algorithm>
#include <iostream>
#include <cmath>
#include <fstream>
#include <gsl/gsl_integration.h>
#include "powerSpectrum.hpp"

using namespace std;

#define USE_GSL
#define TOLERANCE 1e-6
#define NUM_ITERATION 8000

#define POWER_EFSTATHIOU 1
#define HU_WIGGLES 2
#define HU_BARYON 3
#define OLD_POWERSPECTRUM 4
#define POWER_BARDEEN 5
#define POWER_SUGIYAMA 6
#define POWER_BDM 7
#define POWER_TEST 8

#define POWER_SPECTRUM POWER_EFSTATHIOU

namespace Cosmology {

double n = 1.0;
double K0 = 1;
double V0 = 627;

  double CMB_VECTOR[3] = {
    56.759,
    -540.02,
    313.50
  };


// WMAP5
double h = 0.719;
double SIGMA8 = 0.77;
double OMEGA_B = 0.043969;
double OMEGA_C = 0.21259;


// WMAP5-modification
//double h = 0.719;
//double SIGMA8 = 0.77;
//double OMEGA_B = 0;
//double OMEGA_C = 0.21259+0.043969;

// LCDM STRAUSS ?
//double h = 0.67;
//double SIGMA8 = 0.67;
//double OMEGA_B = 0;
//double OMEGA_C = 0.30; 

// SCDM STRAUSS
//double h = 0.5;
//double SIGMA8= 1.05;
//double OMEGA_B = 0;
//double OMEGA_C = 1;

// Sugiyama test
//double h = 0.5;
//double SIGMA8= 0.5;//1.05;
//double OMEGA_B = 0.0125*4;
//double OMEGA_C = 0.1-OMEGA_B;

// HU TEST
//double h = 0.5;
//double SIGMA8 = 0.5;
//double OMEGA_B = 0.09;
//double OMEGA_C = 0.21;

// HDM STRAUSS
//double h = 0.5;
//double SIGMA8 = 0.86;
//double OMEGA_B = 0;
//double OMEGA_C = 1;

// FOR "BEST FIT"
//double h = 0.82;
//double SIGMA8 = 0.76;
//double OMEGA_B = 0.043969;
//double OMEGA_C = 0.15259;


// FOR JUSZKIEWICZ CHECKING (CDM) ! WARNING ! He smoothes 
// with a gaussian filter the density field, i.e. one has
// to multiply P(k) by exp(-k^2 R^2) with R the radius
// of the filter. Dammit !
//double h = 0.5;
//double SIGMA8=1/2.5;
//double OMEGA_B=0.;
//double OMEGA_C=1;
//#define JUSZKIEWICZ_PATCH
//#define RJUSZ 6.0

// (BDM)
//double h = 0.5;
//double SIGMA8=1;
//double OMEGA_B=0.0;
//double OMEGA_C=0.4;


// FOR HU CHECKING
//double h = 0.5;
//double SIGMA8= 1;
//double OMEGA_B=0.09;
//double OMEGA_C=0.21;


double OMEGA_0 = OMEGA_B+OMEGA_C;
double Omega = OMEGA_0;
double Theta_27 = 2.728 / 2.7;
double beta = pow(OMEGA_0, 5./9);
double OmegaEff = OMEGA_0 * h * h;
double Gamma0 = OMEGA_0 * h * h;


/*
 * This is \hat{tophat}
 */
double tophatFilter(double u)
{
  if (u != 0)    
    return 3 / (u*u*u) * (sin(u) - u * cos(u));
  else
    return 1;
}

/*
 * This is \tilde{w}
 */
double externalFilter(double u)
{
  if (u != 0)
    return 1 - sin(u)/u;
  return 0.;
}

double powC(double q, double alpha_c)
{
  return 14.2 / alpha_c + 386 / (1 + 69.9 * pow(q, 1.08));
}

double T_tilde_0(double q, double alpha_c, double beta_c)
{
  double a = log(M_E + 1.8 * beta_c * q);
  return a / ( a  + powC(q, alpha_c) * q * q);
}

double j_0(double x)
{
  if (x == 0)
    return 1.0;
  return sin(x)/x;
}

double powG(double y)
{
  double a = sqrt(1 + y);

  return y * (-6 * a + (2 + 3 * y) *log((a + 1)/(a - 1)));
}

/*
 * This function returns the power spectrum evaluated at k (in Mpc, not in Mpc/h).
 */
double powerSpectrum(double k, double normPower)
{
#if POWER_SPECTRUM == POWER_EFSTATHIOU
  double a = 6.4/Gamma0;
  double b = 3/Gamma0;
  double c = 1.7/Gamma0;
  double nu = 1.13;
  
  double f = (a*k) + pow(b*k,1.5) + pow(c*k,2);

  // EFSTATHIOU  ET AL. (1992)
  return normPower * pow(k,n) * pow(1+pow(f,nu),(-2/nu));
#endif

  // EISENSTEIN ET HU (1998)
  // FULL POWER SPECTRUM WITH BARYONS AND WIGGLES

#if POWER_SPECTRUM == HU_WIGGLES
  // EISENSTEIN ET HU (1998)
  // FULL POWER SPECTRUM WITH BARYONS AND WIGGLES

  double k_silk = 1.6 * pow(OMEGA_B * h * h, 0.52) * pow(OmegaEff, 0.73) * (1 + pow(10.4 * OmegaEff, -0.95));
  double z_eq = 2.50e4 * OmegaEff * pow(Theta_27, -4);
  double s = 44.5 * log(9.83 / OmegaEff) / (sqrt(1 + 10 * pow(OMEGA_B * h * h, 0.75)));
  double f = 1 / (1 + pow(k * s / 5.4, 4));
  double k_eq = 7.46e-2 * OmegaEff * pow(Theta_27, -2);
  double a1 = pow(46.9 * OmegaEff, 0.670) * (1 + pow(32.1 * OmegaEff, -0.532));
  double a2 = pow(12.0 * OmegaEff, 0.424) * (1 + pow(45.0 * OmegaEff, -0.582));
  double alpha_c = pow(a1, -OMEGA_B/ OMEGA_0) * pow(a2, -pow(OMEGA_B / OMEGA_0, 3));  

  double q = k / (13.41 * k_eq);
  double b1_betac = 0.944 * 1/(1 + pow(458 * OmegaEff, -0.708));
  double b2_betac = pow(0.395 * OmegaEff, -0.0266);
  double beta_c = 1/ ( 1 + b1_betac * (pow(OMEGA_C / OMEGA_0, b2_betac) - 1)   );
  double T_c = f * T_tilde_0(q, 1, beta_c) + (1 - f) * T_tilde_0(q, alpha_c, beta_c);

  double b1_zd = 0.313 * pow(OmegaEff, -0.419) * (1 + 0.607 * pow(OmegaEff, 0.674));
  double b2_zd = 0.238 * pow(OmegaEff, 0.223);
  double z_d = 1291 * pow(OmegaEff, 0.251) / (1 + 0.659 * pow(OmegaEff, 0.828)) * (1 + b1_zd * pow(OmegaEff, b2_zd));
  double R_d = 31.5 * OMEGA_B * h * h * pow(Theta_27, -4) * 1e3 / z_d;

  double alpha_b = 2.07 * k_eq * s * pow(1 + R_d, -0.75) * powG((1 + z_eq)/(1 + z_d));
  double beta_b = 0.5 + OMEGA_B / OMEGA_0 + (3 - 2 * OMEGA_B / OMEGA_0) * sqrt(pow(17.2 * OmegaEff, 2) + 1);
  double beta_node = 8.41 * pow(OmegaEff, 0.435);
  double s_tilde = s * pow(1 + pow(beta_node / (k * s), 3), -1./3);

  double T_b = (T_tilde_0(q, 1, 1) / (1 + pow(k * s / 5.2, 2)) + alpha_b / (1 + pow(beta_b / (k * s), 3)) * exp(-pow(k/k_silk, 1.4))) * j_0(k * s_tilde);  

  double T_k = OMEGA_B/OMEGA_0 * T_b + OMEGA_C/OMEGA_0 * T_c;  

  return normPower * pow(k,n) * T_k * T_k;
#endif

  // EISENSTEIN ET AL. (2008), SHAPED POWER SPECTRUM WITH BARYON, WITHOUT WIGGLES
#if POWER_SPECTRUM == HU_BARYON
  double s = 44.5 * log(9.83 / OmegaEff) / (sqrt(1 + 10 * pow(OMEGA_B * h * h, 0.75)));
  double alpha_Gamma = 1 - 0.328 * log(431 * OmegaEff) * OMEGA_B / OMEGA_0 + 0.38 * log(22.3 * OmegaEff) * pow(OMEGA_B / OMEGA_0, 2);
  double GammaEff = OMEGA_0 * h * (alpha_Gamma + (1 - alpha_Gamma)/(1 + pow(0.43 * k * s, 4)));
  double q = k/(h*GammaEff) * pow(Theta_27, 2);
  double L_0 = log(2 * M_E + 1.8 * q);
  double C_0 = 14.2 + 731 / (1 + 62.5 * q);
  double T0 = L_0 / (L_0 + C_0 * q * q);
  
  return normPower * pow(k,n) * T0 * T0;
#endif

#if POWER_SPECTRUM == OLD_POWERSPECTRUM
  // OLD FUNCTION: 
  static const double l = 1/(Omega * h*h);
  static const double alpha = 1.7 * l, beta = 9.0 * pow(l, 1.5), gamma = l*l;
  return normPower * pow(k,n) * pow(1 + alpha * k +  beta * pow(k,1.5) + gamma *k*k,-2);
#endif

#if POWER_SPECTRUM == POWER_SUGIYAMA
  double q = k * Theta_27*Theta_27 / (OmegaEff * exp(-OMEGA_B - sqrt(h/0.5)*OMEGA_B/OMEGA_0));
  double L0 = log(2*M_E + 1.8 * q);
  double C0 = 14.2 + 731 / (1 + 62.5 * q);
  double T_k = L0 / (L0 + C0 * q*q);

      //  double poly = 1 + 3.89 * q + pow(16.1*q,2) + pow(5.46*q,3) + pow(6.71*q,4);
      //  double T_k = log(1+2.34*q)/(2.34*q) * pow(poly,-0.25);

  return normPower * pow(k,n) * T_k * T_k;

#endif

#if POWER_SPECTRUM == POWER_BARDEEN
  double q = k / (OmegaEff);

  double poly = 1 + 3.89 * q + pow(16.1*q,2) + pow(5.46*q,3) + pow(6.71*q,4);
  double T_k = log(1+2.34*q)/(2.34*q) * pow(poly,-0.25);

  return normPower * pow(k,n) * T_k * T_k;

#endif

#if POWER_SPECTRUM == POWER_BDM

  k /= h*h;
  double alpha1 = 190;
  double Gmu = 4.6;
  double alpha2 = 11.5;
  double alpha3 = 11;
  double alpha4 = 12.55;
  double alpha5 = 0.0004;
  return normPower*k*alpha1*alpha1*Gmu*Gmu/(1+(alpha2*k)+pow(alpha3*k,2)+pow(alpha4*k,3))*pow(1+pow(alpha5/k,2), -2);
#endif

#if POWER_SPECTRUM == POWER_TEST
  return 1/(1+k*k);
#endif
}

/*
 * This function computes the normalization of the power spectrum. It requests
 * a sigma8 (density fluctuations within 8 Mpc/h)
 */
double gslPowSpecNorm(double k, void *params)
{
  double f = tophatFilter(k*8.0/h);
  return powerSpectrum(k, 1.0)*k*k*f*f;
}

double computePowSpecNorm(double sigma8)
{
  int Nsteps = 30000;
  double normVal = 0;

#ifndef USE_GSL
  for (int i = 1; i <= Nsteps; i++)
    {
      double t = i * 1.0/(Nsteps+1);
      // Change of variable !
      double k = (1-t)/t * K0;

      // The filter
      double filter_val = tophatFilter(k*8.0/h);
      // The powerspectrum
      double powVal = powerSpectrum(k, 1.0);
      // Multiply by the tophat filter
      powVal *= filter_val*filter_val;
     
      powVal *= k*k;
      
      // Account for change of variable
      powVal /= (t*t);
  
      // Integrate !
      normVal += powVal;
    }

  normVal /= 2*M_PI*M_PI;

  // The dt element
  normVal *= 1.0/(Nsteps+1) * K0;
#else
  double abserr;
  gsl_integration_workspace *w = gsl_integration_workspace_alloc(NUM_ITERATION);
  gsl_function f;

  f.function = gslPowSpecNorm;
  f.params = 0;

  gsl_integration_qagiu(&f, 0, 0, TOLERANCE, NUM_ITERATION, w, &normVal, &abserr);

  gsl_integration_workspace_free(w);

  normVal /= (2*M_PI*M_PI);
#endif

  return sigma8*sigma8/normVal;
}

/*
 * This function computes the variance of the Local Group velocity components
 * for a survey which depth is topHatRad1 (in Mpc/h). This variance should
 * be multiplied by (H \beta)^2 to be equal to a velocity^2.
 */
double gslVariance(double k, void *params)
{
  double R1 = *(double *)params;
  double f = externalFilter(k * R1 / h);

  double a = f*f;

#ifdef JUSZKIEWICZ_PATCH
  a *= exp(-k*k*(RJUSZ*RJUSZ/(h*h)));
#endif
  a *= powerSpectrum(k, 1.0);

  return a;
}

double computeVariance(double powNorm,  double topHatRad1)
{
  int Nsteps = 100000;
  double varVal = 0;

#ifndef USE_GSL
  
  for (int i = 1; i <= Nsteps; i++)
    {
      double t = i * 1.0/(Nsteps+1);
      // Change of variable !
      double k = (1-t)/t * K0;
      double powVal = powerSpectrum(k, powNorm);

      double filter1Val = externalFilter(k*topHatRad1/h);

#ifdef JUSZKIEWICZ_PATCH
      powVal *= exp(-k*k*(RJUSZ*RJUSZ/(h*h)));
#endif
      powVal *= filter1Val*filter1Val;
      
      powVal /= (t*t);

      varVal += powVal;
    }

  varVal *= 1.0/(Nsteps) * K0;
  varVal *= 1.0/(6*M_PI*M_PI);

#else

  double abserr;
  gsl_integration_workspace *w = gsl_integration_workspace_alloc(NUM_ITERATION);
  gsl_function f;
  
  f.function = gslVariance;
  f.params = &topHatRad1;

  gsl_integration_qagiu(&f, 0, 0, TOLERANCE, NUM_ITERATION, w, &varVal, &abserr);

  gsl_integration_workspace_free(w);

  varVal *= powNorm/(6*M_PI*M_PI);

#endif

  return varVal;
}

/*
 * This function computes the same quantity as computeVariance but
 * for a survey infinitely deep.
 */
double gslVarianceZero(double k, void *params)
{
  double a = 1.0;

#ifdef JUSZKIEWICZ_PATCH
  a *= exp(-k*k*(RJUSZ*RJUSZ/(h*h)));
#endif
  a *= powerSpectrum(k, 1.0);

  return a;
}

double computeVarianceZero(double powNorm)
{
  int Nsteps = 100000;
  double varVal = 0;
  
#ifndef USE_GSL

  for (int i = 1; i <= Nsteps; i++)
    {
      double t = i * 1.0/(Nsteps+1);
      // Change of variable !
      double k = (1-t)/t * K0;
      double powVal = powerSpectrum(k, powNorm);
      
#ifdef JUSZKIEWICZ_PATCH
      powVal *= exp(-k*k*(RJUSZ*RJUSZ/h*h));
#endif

      powVal /= (t*t);

      varVal += powVal;
    }

  varVal *= 1.0/(Nsteps+1) * K0;
  varVal *= 1.0/(6*M_PI*M_PI);

#else
 
  double abserr;
  gsl_integration_workspace *w = gsl_integration_workspace_alloc(NUM_ITERATION);
  gsl_function f;
  
  f.function = gslVarianceZero;
  f.params = 0;

  gsl_integration_qagiu(&f, 0, 0, TOLERANCE, NUM_ITERATION, w, &varVal, &abserr);

  gsl_integration_workspace_free(w);

  varVal *= powNorm/(6*M_PI*M_PI);
#endif

  return varVal;
}

/*
 * This function computes the correlation between the infinitely deep
 * velocity of the Local Group and the one estimated from a survey
 * which depth is topHatRad1.
 * This corresponds to \gamma.
 * This quantity must be multiplied by H \beta to be equal to a velocity^2.
 */
double gslCorrel(double k, void *params)
{
  double R1 = ((double *)params)[0];
  double a = externalFilter(k * R1 / h);// * externalFilter(k * R2 / h);

#ifdef JUSZKIEWICZ_PATCH
  a *= exp(-k*k*(RJUSZ*RJUSZ/(h*h)));
#endif
  a *= powerSpectrum(k, 1.0);

  return a;
}

double gslCorrelBis(double t, void *params)
{
  double k = (1-t)/t;
  double v = gslCorrel(k, params);

  return v/(t*t);
}

double gslCorrel2(double k, void *params)
{
  double R1 = ((double *)params)[0];
  double R2 = ((double *)params)[1];
  double a = externalFilter(k * R1 / h) * externalFilter(k * R2 / h);        
  
#ifdef JUSZKIEWICZ_PATCH
  a *= exp(-k*k*(RJUSZ*RJUSZ/(h*h)));
#endif
  a *= powerSpectrum(k, 1.0) ;
  
  return a;
}

double gslCorrel2bis(double t, void *params)
{
  double k = (1-t)/t;
  double v = gslCorrel2(k, params);

  return v/(t*t);
}


double computeCorrel(double powNorm,  double topHatRad1)
{
  int Nsteps = 100000;
  double varVal = 0;
  
#ifndef USE_GSL

  for (int i = 1; i <= Nsteps; i++)
    {
      double t = i * 1.0/(Nsteps+1);
      // Change of variable !
      double k = (1-t)/t * K0;
      double powVal = powerSpectrum(k, powNorm);

      double filter1Val = externalFilter(k*topHatRad1/h);

#ifdef JUSZKIEWICZ_PATCH
      powVal*=exp(-k*k*(RJUSZ*RJUSZ/h*h));
#endif

      powVal *= filter1Val;
      
      powVal /= (t*t);

      varVal += powVal;
    }

  varVal *= 1.0/(Nsteps) * K0;
  varVal *= 1.0/(6*M_PI*M_PI);
#else

  double abserr;
  gsl_integration_workspace *w = gsl_integration_workspace_alloc(NUM_ITERATION);
  gsl_function f;
  
  f.params = &topHatRad1;

#if 1
  f.function = gslCorrelBis;
  gsl_integration_qag (&f, 0, 1, 0, TOLERANCE, NUM_ITERATION, GSL_INTEG_GAUSS61, w, &varVal, &abserr);
#else
  f.function = gslCorrel;
  gsl_integration_qagiu(&f, 0, 0, TOLERANCE, NUM_ITERATION, w, &varVal, &abserr);
#endif

  gsl_integration_workspace_free(w);

  varVal *= powNorm/(6*M_PI*M_PI);

#endif
  return varVal;
}

/*
 * This function computes the correlation between the infinitely deep
 * velocity of the Local Group and the one estimated from a survey
 * which depth is topHatRad1.
 * This corresponds to \gamma.
 * This quantity must be multiplied by H \beta to be equal to a velocity^2.
 */
double computeCorrel2(double powNorm,  double topHatRad1, double topHatRad2)
{
  int Nsteps = 100000;
  double varVal = 0;
#ifndef USE_GSL
  
  for (int i = 1; i <= Nsteps; i++)
    {
      double t = i * 1.0/(Nsteps+1);
      // Change of variable !
      double k = (1-t)/t * K0;
      double powVal = powerSpectrum(k, powNorm);

      double filter1Val = externalFilter(k*topHatRad1/h);
      double filter2Val = externalFilter(k*topHatRad2/h);

      powVal *= filter1Val * filter2Val;
      
      powVal /= (t*t);

      varVal += powVal;
    }

  varVal *= 1.0/(Nsteps) * K0;
  varVal *= 1.0/(6*M_PI*M_PI);

#else

  double abserr;
  gsl_integration_workspace *w = gsl_integration_workspace_alloc(NUM_ITERATION);
  gsl_function f;
  double rads[] = {topHatRad1, topHatRad2 };

  f.params = rads;

#if 1  
  f.function = gslCorrel2bis;
  gsl_integration_qag (&f, 0, 1, 0, TOLERANCE, NUM_ITERATION, GSL_INTEG_GAUSS61, w, &varVal, &abserr);
#else
  f.function = gslCorrel2;
  gsl_integration_qagiu(&f, 0, 0, TOLERANCE, NUM_ITERATION, w, &varVal, &abserr);
#endif

  gsl_integration_workspace_free(w);

  varVal *= powNorm/(6*M_PI*M_PI);

#endif

  return varVal;
}


  void updateCosmology()
  {
    OMEGA_0 = OMEGA_B+OMEGA_C;
    Omega = OMEGA_0;
    Theta_27 = 2.728 / 2.7;
    beta = pow(OMEGA_0, 5./9);
    OmegaEff = OMEGA_0 * h * h;
    Gamma0 = OMEGA_0 * h * h;

#if 0
    cout << "Cosmology is :" << endl
	 << "  O0=" << OMEGA_0 << " Theta=" << Theta_27 << " beta=" << beta << " h=" << h << " G0=" << Gamma0 << endl
	 << "  OmegaB=" << OMEGA_B << " Omega_C=" << OMEGA_C << endl;
#endif
  }

  double vvCorrection(double P_deltadelta, double k)
  {
    static const double alpha0 = -12480.5, alpha1 = 1.824, alpha2 = 2165.87, alpha3=1.796;
    if (k > 0.3)
      return 0;
    double r =(alpha0*sqrt(P_deltadelta) + alpha1*P_deltadelta*P_deltadelta)/(alpha2 + alpha3*P_deltadelta);
    assert(P_deltadelta > 0);

    if (r < 0)
      return 0;
    return r;
  }

};
