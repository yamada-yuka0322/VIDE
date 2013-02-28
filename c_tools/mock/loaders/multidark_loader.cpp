#include <cassert>
#include <string>
#include <fstream>
#include <iostream>
#include <CosmoTool/loadRamses.hpp>
#include <CosmoTool/fortran.hpp>
#include "simulation_loader.hpp"

using namespace std;
using namespace CosmoTool;

class MultiDarkLoader: public SimulationLoader
{
protected:
  SimuData *header;
  string darkname;
  SimulationPreprocessor *preproc;
public:
  MultiDarkLoader(const std::string& name, SimuData *h, SimulationPreprocessor *p)
    : preproc(p), darkname(name), header(h)
  {
  }

  ~MultiDarkLoader()
  {
    delete header;
  }

  int num_files()
  {
    return 1;
  }

  SimuData *getHeader()
  {
    return header;
  }

  SimuData *loadFile(int id)
  {
    if (id != 0)
      return 0;

    ifstream fp(darkname.c_str());
    SimuData *simu = new SimuData;

    fp >> simu->BoxSize >> simu->Omega_M >> simu->Hubble >> simu->time >> simu->NumPart;
    simu->time = 1./(1.+simu->time); // convert to scale factor
    simu->TotalNumPart = simu->NumPart;
    simu->Omega_Lambda = 1.0 - simu->Omega_M;

    long estimated = (preproc == 0) ? simu->NumPart : preproc->getEstimatedPostprocessed(simu->NumPart);
    long allocated = estimated; 

    for (int k = 0; k < 3; k++)
      simu->Pos[k] = new float[allocated];
    simu->Vel[2] = new float[allocated];
    simu->Id = new long[allocated];
    long *uniqueID = new long[allocated];
    long *index = new long[allocated];

    double tempData;

    simu->new_attribute("uniqueID", uniqueID, delete_adaptor<long>);
    simu->new_attribute("index", index, delete_adaptor<long>);

    cout << "loading multidark particles" << endl;
    long actualNumPart = 0;

    for (long i = 0; i < simu->NumPart; i++) {
      SingleParticle p;

      fp >> p.ID >> p.Pos[0] >> p.Pos[1]
         >> p.Pos[2] >> p.Vel[2] >> tempData >> tempData;

      if (p.ID == -99 && 
          p.Pos[0] == -99 && p.Pos[1] == -99 && 
          p.Pos[2] == -99 && p.Vel[2] == -99) {
        break;

      if (preproc != 0 && !preproc->accept(p))
        continue;

      copyParticleToSimu(p, simu, actualNumPart);
      uniqueID[actualNumPart]= p.ID;
      index[actualNumPart] = i;
      actualNumPart++;
      if (actualNumPart == allocated)
        {
          allocated += (estimated+9)/10;
          reallocSimu(simu, allocated); 
          reallocArray(uniqueID, allocated, actualNumPart);
          reallocArray(index, allocated, actualNumPart);
        }
      }
    }
    applyTransformations(simu);
    simu->NumPart = actualNumPart;
    simu->TotalNumPart = actualNumPart;
    return simu;
  }
};

SimulationLoader *multidarkLoader(const string& multidarkname, SimulationPreprocessor *p)
{
  SimuData *header;
  int actualNumPart;
  ifstream fp(multidarkname.c_str());

  cout << "opening multidark file " << multidarkname << endl;
  if (!fp)
    {
      cout << "could not open file!" << endl;
      return 0;
    }

  header = new SimuData();
  fp >> header->BoxSize >> header->Omega_M >> header->Hubble >> header->time >> header->NumPart;
  
  header->time = 1./(1.+header->time); // convert to scale factor
  header->TotalNumPart = header->NumPart;
  header->Omega_Lambda = 1.0 - header->Omega_M;

  return new MultiDarkLoader(multidarkname, header, p);
}

  

