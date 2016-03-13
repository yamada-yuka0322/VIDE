/*+
This is CosmoTool (./src/load_data.hpp) -- Copyright (C) Guilhem Lavaux (2007-2014)

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
    int flag_doubleprecision;	
    int flag_ic_info;            
    float lpt_scalingfactor;
    char fill[18];		/*!< fills to 256 Bytes */
    char names[15][2];
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
