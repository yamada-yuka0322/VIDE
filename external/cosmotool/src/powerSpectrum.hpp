/*+
This is CosmoTool (./src/powerSpectrum.hpp) -- Copyright (C) Guilhem Lavaux (2007-2014)

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
