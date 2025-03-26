/*+
This is CosmoTool (./src/cosmopower.hpp) -- Copyright (C) Guilhem Lavaux (2007-2014)

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

#ifndef _COSMOPOWER_HPP
#define _COSMOPOWER_HPP

namespace CosmoTool {

  struct TF_Transfer;

  class CosmoPower
  {
  public:
    // PRIMARY VARIABLES
    double n;
    double K0;
    double V_LG_CMB;

    double CMB_VECTOR[3];

    double h;
    double SIGMA8;
    double OMEGA_B;
    double OMEGA_C;
    double omega_B;
    double omega_C;
    double Theta_27;

    // DERIVED VARIABLES
    double OMEGA_0;
    double Omega;
    double beta;
    double OmegaEff;
    double Gamma0;
    double normPower;

    struct EHuParams {
      double k_silk;
      double s;
      double k_eq;
      double alpha_b, beta_b;
      double alpha_c, beta_c;
      double beta_node;
    };

    EHuParams ehu;
    TF_Transfer *ehu_params;

    enum CosmoFunction {
      POWER_EFSTATHIOU,
      HU_WIGGLES,
      PRIMORDIAL_PS,
      MATTER_TK,
      HU_BARYON,
      OLD_POWERSPECTRUM,
      POWER_BARDEEN,
      POWER_SUGIYAMA,
      POWER_BDM,
      POWER_TEST,
      HU_WIGGLES_ORIGINAL
    };

    CosmoPower();
    ~CosmoPower();

    void setFunction(CosmoFunction f);

    void updateCosmology();
    void updatePhysicalCosmology();
    void normalize(double k_min = -1, double k_max = -1);
    void setNormalization(double A_K);
    void updateHuWigglesConsts();
    void updateHuWigglesOriginal();

    double eval_theta_theta(double k);
    double power(double k);

    double integrandNormalize(double k);
  private:
    double (CosmoPower::*eval)(double);

    double powerEfstathiou(double k);
    double powerHuWiggles(double k);
    double primordialPowerSpectrum(double k);
    double matterTransferFunctionHu(double k);
    double powerHuBaryons(double k);
    double powerOld(double k);
    double powerBardeen(double k);
    double powerSugiyama(double k);
    double powerBDM(double k);
    double powerTest(double k);
    double powerHuWigglesOriginal(double k);
  };

};

#endif
