#ifndef _LOAD_GADGET_DATA_HPP
#define _LOAD_GADGET_DATA_HPP

#include "config.hpp"

namespace CosmoTool {

  struct GadgetHeader
  {
    int      npart[6];
    double   mass[6];
    double   time;
    double   redshift;
    int      flag_sfr;
    int      flag_feedback;
    int      npartTotal[6];
    int      flag_cooling;
    int      num_files;
    double   BoxSize;
    double   Omega0;
    double   OmegaLambda;
    double   HubbleParam; 
    char     fill[256- 6*4- 6*8- 2*8- 2*4- 6*4- 2*4 - 4*8];  /* fills to 256 Bytes */
  };

  struct ParticleState
  {
    float  Pos[3];
    float  Init[3];
    float  Vel[3];
    float  VelInit[3];
    float  Mass;
    int    Type;
    
    float  Rho, U, Temp, Ne;
    int Id;
  };
  
  struct GadgetData {
    GadgetHeader header;
    ParticleState *particles;
    int NumPart;
    int ntot_withmasses;
  };
  
  struct ParticleSetHeader {
    int Ntypes;
    float BoxSize;
    float RedShift;
    char header[256 - 4 - 2*4];
  };

  struct ParticleSet {
    ParticleSetHeader header;
    // Particle description
    int *Npart;
    ParticleState **particles;
  };

  struct PurePositionData {
    unsigned int NumPart;
    double BoxSize;
    double hubble;
    FCoordinates *pos;
  };

  struct PhaseSpaceData {
    unsigned int NumPart;
    double hubble;
    double BoxSize;
    FCoordinates *pos;
    FCoordinates *vel;
  };

  struct PhaseSpaceDataID {
    unsigned int NumPart;
    double hubble;
    double BoxSize;
    FCoordinates *pos;
    FCoordinates *vel;
    int *ID;
  };

  void writeGadget(GadgetData *data, const char *fname);
  GadgetData *loadGadget(const char *fname);
  void freeGadget(GadgetData *data);

  GadgetData *loadSimulationData(const char *fname);
  GadgetData *loadHydra(const char *fname);

  void writePersoSet(ParticleSet *set, const char *fname);
  ParticleSet *loadPersoSet(const char *fname);
  void freePersoSet(ParticleSet *set);

};

#endif
