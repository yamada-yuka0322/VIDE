/*+
This is CosmoTool (./src/cosmopower.cpp) -- Copyright (C) Guilhem Lavaux (2007-2014)

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
#include "cosmopower.hpp"
#include "tf_fit.hpp"

using namespace std;
using namespace CosmoTool;

#define USE_GSL
#define TOLERANCE 1e-6
#define NUM_ITERATION 8000


CosmoPower::CosmoPower()
{
  eval = &CosmoPower::powerEfstathiou;
  
  n = 1.0;
  K0 = 1;
  V_LG_CMB = 627;

  CMB_VECTOR[0] = 56.759;
  CMB_VECTOR[1] = -540.02;
  CMB_VECTOR[2] = 313.50;

  h = 0.719;
  SIGMA8 = 0.77;
  OMEGA_B = 0.043969;
  OMEGA_C = 0.21259;

  Theta_27 = 2.728/2.7;

  ehu_params = 0;
  
  updateCosmology();
}

CosmoPower::~CosmoPower()
{
  if (ehu_params)
    delete ehu_params;
}

/*
 * This is \hat{tophat}
 */
static double tophatFilter(double u)
{
  if (u != 0)    
    return 3 / (u*u*u) * (sin(u) - u * cos(u));
  else
    return 1;
}

static double powC(double q, double alpha_c)
{
  return 14.2 / alpha_c + 386 / (1 + 69.9 * pow(q, 1.08));
}

static double T_tilde_0(double q, double alpha_c, double beta_c)
{
  static const double c_E = M_E;
  double a = log(c_E + 1.8 * beta_c * q);
  return a / ( a  + powC(q, alpha_c) * q * q);
}

static double j_0(double x)
{
  if (x == 0)
    return 1.0;
  return sin(x)/x;
}

static double powG(double y)
{
  double a = sqrt(1 + y);

  return y * (-6 * a + (2 + 3 * y) *log((a + 1)/(a - 1)));
}
 
double CosmoPower::powerEfstathiou(double k)
{
  double a = 6.4/Gamma0;
  double b = 3/Gamma0;
  double c = 1.7/Gamma0;
  double nu = 1.13;
  
  double f = (a*k) + pow(b*k,1.5) + pow(c*k,2);

  // EFSTATHIOU  ET AL. (1992)
  return normPower * pow(k,n) * pow(1+pow(f,nu),(-2/nu));
}

void CosmoPower::updateHuWigglesOriginal()
{
  if (ehu_params == 0)
    ehu_params = new TF_Transfer();
  
  ehu_params->TFset_parameters( (OMEGA_C+OMEGA_B)*h*h,
                     OMEGA_B/(OMEGA_C+OMEGA_B), Theta_27*2.7);
}

void CosmoPower::updateHuWigglesConsts()
{
  double f_b = OMEGA_B / OMEGA_0;
  double f_c = OMEGA_C / OMEGA_0;

  double k_silk = 1.6 * pow(OMEGA_B * h * h, 0.52) * pow(OmegaEff, 0.73) * (1 + pow(10.4 * OmegaEff, -0.95));
  double z_eq = 2.50e4 * OmegaEff * pow(Theta_27, -4);
  //double s = 44.5 * log(9.83 / OmegaEff) / (sqrt(1 + 10 * pow(OMEGA_B * h * h, 0.75)));
  double k_eq = 7.46e-2 * OmegaEff * pow(Theta_27, -2);

  double b1_zd = 0.313 * pow(OmegaEff, -0.419) * (1 + 0.607 * pow(OmegaEff, 0.674));
  double b2_zd = 0.238 * pow(OmegaEff, 0.223);
  double z_d = 1291 * pow(OmegaEff, 0.251) / (1 + 0.659 * pow(OmegaEff, 0.828)) * (1 + b1_zd * pow(OMEGA_B*h*h, b2_zd));

  double R_d = 31.5 * OMEGA_B * h * h * pow(Theta_27, -4) * 1e3 / z_d;
  double Req = 31.5 * OMEGA_B * h * h * pow(Theta_27, -4) * 1e3 / z_eq;

  double s = 2./(3.*k_eq) * sqrt(6/Req) * log((sqrt(1 + R_d) + sqrt(R_d + Req))/(1 + sqrt(Req)));

  double a1 = pow(46.9 * OmegaEff, 0.670) * (1 + pow(32.1 * OmegaEff, -0.532));
  double a2 = pow(12.0 * OmegaEff, 0.424) * (1 + pow(45.0 * OmegaEff, -0.582));
  double alpha_c = pow(a1, -f_b) * pow(a2, -pow(f_b, 3));  

  double b1_betac = 0.944 * 1/(1 + pow(458 * OmegaEff, -0.708));
  double b2_betac = pow(0.395 * OmegaEff, -0.0266);
  double beta_c = 1/ ( 1 + b1_betac * (pow(f_c, b2_betac) - 1)   );

  double alpha_b = 2.07 * k_eq * s * pow(1 + R_d, -0.75) * powG((1 + z_eq)/(1 + z_d));
  double beta_b = 0.5 + f_b + (3 - 2 * f_b) * sqrt(pow(17.2 * OmegaEff, 2) + 1);
  double beta_node = 8.41 * pow(OmegaEff, 0.435);
  
  ehu.k_silk = k_silk;
  ehu.s = s;
  ehu.k_eq = k_eq;
  ehu.alpha_c = alpha_c;
  ehu.beta_c = beta_c;
  
  ehu.alpha_b = alpha_b;
  ehu.beta_b = beta_b;
  ehu.beta_node = beta_node;
}

double CosmoPower::powerHuWiggles(double k)
{
  // EISENSTEIN ET HU (1998)
  // FULL POWER SPECTRUM WITH BARYONS AND WIGGLES

  double k_silk = ehu.k_silk;
  double s = ehu.s;
  double k_eq = ehu.k_eq;
  double alpha_c = ehu.alpha_c;
  double beta_c = ehu.beta_c;
  double alpha_b = ehu.alpha_b;
  double beta_b = ehu.beta_b;
  double xx = k * s;  

  double s_tilde = s * pow(1 + pow(ehu.beta_node / (xx), 3), -1./3);
  
  double f = 1 / (1 + pow(xx / 5.4, 4));
  double q = k / (13.41 * k_eq);
  double T_c = f * T_tilde_0(q, 1, beta_c) + (1 - f) * T_tilde_0(q, alpha_c, beta_c);

  double T_b = ( 
     T_tilde_0(q, 1, 1) / (1 + pow(xx / 5.2, 2)) + 
     alpha_b / (1 + pow(beta_b / xx, 3)) * exp(-pow(k/k_silk, 1.4))
     ) * j_0(k * s_tilde);  

  double T_k = OMEGA_B/OMEGA_0 * T_b + OMEGA_C/OMEGA_0 * T_c;  

  return normPower * pow(k,n) * T_k * T_k;
}

double CosmoPower::powerHuBaryons(double k)
{
  double s = 44.5 * log(9.83 / OmegaEff) / (sqrt(1 + 10 * pow(OMEGA_B * h * h, 0.75)));
  double alpha_Gamma = 1 - 0.328 * log(431 * OmegaEff) * OMEGA_B / OMEGA_0 + 0.38 * log(22.3 * OmegaEff) * pow(OMEGA_B / OMEGA_0, 2);
  double GammaEff = OMEGA_0 * h * (alpha_Gamma + (1 - alpha_Gamma)/(1 + pow(0.43 * k * s, 4)));
  double q = k/(h*GammaEff) * pow(Theta_27, 2);
  double L_0 = log(2 * M_E + 1.8 * q);
  double C_0 = 14.2 + 731 / (1 + 62.5 * q);
  double T0 = L_0 / (L_0 + C_0 * q * q);
  
  return normPower * pow(k,n) * T0 * T0;
}

double CosmoPower::powerOld(double k)
{
  static const double l = 1/(Omega * h*h);
  static const double alpha = 1.7 * l, beta = 9.0 * pow(l, 1.5), gamma = l*l;
  return normPower * pow(k,n) * pow(1 + alpha * k +  beta * pow(k,1.5) + gamma *k*k,-2);
}

double CosmoPower::powerSugiyama(double k)
{
  double q = k * Theta_27*Theta_27 / (OmegaEff * exp(-OMEGA_B - sqrt(h/0.5)*OMEGA_B/OMEGA_0));
  double L0 = log(2*M_E + 1.8 * q);
  double C0 = 14.2 + 731 / (1 + 62.5 * q);
  double T_k = L0 / (L0 + C0 * q*q);

  return normPower * pow(k,n) * T_k * T_k;
}

double CosmoPower::powerBardeen(double k)
{
  double q = k / (OmegaEff);

  double poly = 1 + 3.89 * q + pow(16.1*q,2) + pow(5.46*q,3) + pow(6.71*q,4);
  double T_k = log(1+2.34*q)/(2.34*q) * pow(poly,-0.25);

  return normPower * pow(k,n) * T_k * T_k;
}

double CosmoPower::powerBDM(double k)
{
  k /= h*h;
  double alpha1 = 190;
  double Gmu = 4.6;
  double alpha2 = 11.5;
  double alpha3 = 11;
  double alpha4 = 12.55;
  double alpha5 = 0.0004;
  return normPower*k*alpha1*alpha1*Gmu*Gmu/(1+(alpha2*k)+pow(alpha3*k,2)+pow(alpha4*k,3))*pow(1+pow(alpha5/k,2), -2);
}

double CosmoPower::powerTest(double k)
{
  return normPower;//1/(1+k*k);
}

/*
 * This function computes the normalization of the power spectrum. It requests
 * a sigma8 (density fluctuations within 8 Mpc/h)
 */
static double gslPowSpecNorm(double k, void *params)
{
  CosmoPower *c = (CosmoPower *)params;

  return c->integrandNormalize(k);
}

double CosmoPower::integrandNormalize(double x)
{
  double k = (1-x)/x;
  double f = tophatFilter(k*8.0/h);
  return power(k)*k*k*f*f/(x*x);
}

void CosmoPower::normalize(double k_min, double k_max)
{
  double normVal = 0;
  double abserr;
  gsl_integration_workspace *w = gsl_integration_workspace_alloc(NUM_ITERATION);
  gsl_function f;
  double x_min = 0, x_max = 1;

  if (k_max > 0)
    x_min = 1/(1+k_max);
  if (k_min > 0)
    x_max = 1/(1+k_min);

  f.function = gslPowSpecNorm;
  f.params = this;

  normPower = 1;

  ofstream ff("PP_k.txt");
  for (int i = 0; i < 100; i++)
    {
       double k = pow(10.0, 8.0*i/100.-4);
       ff << k << " " << power(k) << endl;
    }

 // gsl_integration_qagiu(&f, 0, 0, TOLERANCE, NUM_ITERATION, w, &normVal, &abserr);
  gsl_integration_qag(&f, x_min, x_max, 0, TOLERANCE, NUM_ITERATION, GSL_INTEG_GAUSS61, w, &normVal, &abserr);
  gsl_integration_workspace_free(w);

  normVal /= (2*M_PI*M_PI);

  normPower =  SIGMA8*SIGMA8/normVal;
}

void CosmoPower::updateCosmology()
{
  OMEGA_0 = OMEGA_B+OMEGA_C;
  Omega = OMEGA_0;
  beta = pow(OMEGA_0, 5./9);
  OmegaEff = OMEGA_0 * h * h;
  Gamma0 = OMEGA_0 * h * h;
  omega_B = OMEGA_B * h * h;
  omega_C = OMEGA_C * h * h;
}

void CosmoPower::updatePhysicalCosmology()
{
  OMEGA_B = omega_B / (h*h);
  OMEGA_C = omega_C / (h*h);
  OMEGA_0 = Gamma0 / (h*h);
  beta = pow(OMEGA_0, 5./9);
}

double CosmoPower::eval_theta_theta(double k)
{
  // Jennings (2012) fit
  double P_deltadelta = power(k);

  static const double alpha0 = -12480.5, alpha1 = 1.824, alpha2 = 2165.87, alpha3=1.796;
  if (k > 0.3)
    return 0;
  double r =(alpha0*sqrt(P_deltadelta) + alpha1*P_deltadelta*P_deltadelta)/(alpha2 + alpha3*P_deltadelta);
  assert(P_deltadelta > 0);
  
  if (r < 0)
    return 0;
  return r;
}

double CosmoPower::power(double k)
{
  return (this->*eval)(k);
}

double CosmoPower::powerHuWigglesOriginal(double k)
{
  float tfb, tfc;
  
  ehu_params->TFfit_onek(k, &tfb, &tfc);

  double T_k = OMEGA_B/OMEGA_0 * tfb + OMEGA_C/OMEGA_0 * tfc;  

  return normPower * pow(k,n) * T_k * T_k;
}

void CosmoPower::setFunction(CosmoFunction f)
{
  switch (f)
    {
    case POWER_EFSTATHIOU:
      eval = &CosmoPower::powerEfstathiou;
      break;
    case HU_WIGGLES:
      updateHuWigglesConsts();
      eval = &CosmoPower::powerHuWiggles;
      break;
    case HU_WIGGLES_ORIGINAL:
      updateHuWigglesOriginal();
      eval = &CosmoPower::powerHuWigglesOriginal;
      break;
    case HU_BARYON:
      eval = &CosmoPower::powerHuBaryons;
      break;
    case OLD_POWERSPECTRUM:
      eval = &CosmoPower::powerOld;
      break;
    case POWER_BARDEEN:
      eval = &CosmoPower::powerBardeen;
      break;
    case POWER_SUGIYAMA:
      eval = &CosmoPower::powerSugiyama;
      break;
    case POWER_BDM:
      eval = &CosmoPower::powerBDM;
      break;
    case POWER_TEST:
      eval = &CosmoPower::powerTest;
      break;
    default:
      abort();
    }
}

void CosmoPower::setNormalization(double A_K)
{
  normPower = A_K;///power(0.002);
}

