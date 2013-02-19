#ifndef _POWERSPECTRUM_HPP
#define _POWERSPECTRUM_HPP

namespace Cosmology {

  extern double n;
  extern double K0;
  extern double V0;

  extern double CMB_VECTOR[3];

  // WMAP5
  extern double h;
  extern double SIGMA8;
  extern double OMEGA_B;
  extern double OMEGA_C;

  extern double OMEGA_0;
  extern double Omega;
  extern double Theta_27;
  extern double beta;
  extern double OmegaEff;
  extern double Gamma0;

  void updateCosmology();
  double tophatFilter(double u);
  double powerSpectrum(double k, double normPower);
  double computePowSpecNorm(double sigma8);
  double computeVariance(double powNorm,  double topHatRad1);
  double computeVarianceZero(double powNorm);
  double computeCorrel(double powNorm,  double topHatRad1);
  double computeCorrel2(double powNorm,  double topHatRad1, double topHatRad2);

  double vvCorrection(double P_deltadelta, double k);
};

#endif
