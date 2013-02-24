#ifndef _COSMOPOWER_HPP
#define _COSMOPOWER_HPP

namespace CosmoTool {

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

    enum CosmoFunction
      {
	POWER_EFSTATHIOU,
	HU_WIGGLES,
	HU_BARYON,
	OLD_POWERSPECTRUM,
	POWER_BARDEEN,
	POWER_SUGIYAMA,
	POWER_BDM,
	POWER_TEST
      };

    CosmoPower();

    void setFunction(CosmoFunction f);

    void updateCosmology();
    void updatePhysicalCosmology();
    void normalize();
    void setNormalization(double A_K);

 
    double eval_theta_theta(double k);
    double power(double k);

    double integrandNormalize(double k);
  private:
    double (CosmoPower::*eval)(double);

    double powerEfstathiou(double k);
    double powerHuWiggles(double k);
    double powerHuBaryons(double k);
    double powerOld(double k);
    double powerBardeen(double k);
    double powerSugiyama(double k);
    double powerBDM(double k);
    double powerTest(double k);

  };

};

#endif
