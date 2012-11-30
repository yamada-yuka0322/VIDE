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
public:
  MultiDarkLoader(const std::string& name, SimuData *h)
    : darkname(name), header(h)
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

    for (int k = 0; k < 3; k++)
      simu->Pos[k] = new float[simu->NumPart];
    simu->Vel[2] = new float[simu->NumPart];
    simu->Id = new int[simu->NumPart];
    long *uniqueID = new long[simu->NumPart];

    simu->new_attribute("uniqueID", uniqueID, delete_adaptor<long>);

    cout << "loading multidark particles" << endl;
    long actualNumPart = 0;
    for (int i = 0; i < simu->NumPart; i++) {

      fp >> simu->Id[i] >> simu->Pos[0][i] >> simu->Pos[1][i]
	 >> simu->Pos[2][i] >> simu->Vel[2][i];

      uniqueID[i] = simu->Id[i];

      if (simu->Id[i] == -99 && 
          simu->Pos[0][i] == -99 && simu->Pos[1][i] == -99 && 
          simu->Pos[2][i] == -99 && simu->Vel[2][i] == -99) {
        break;
      } else {
        actualNumPart++;
      }
    }

    simu->NumPart = actualNumPart;
    simu->TotalNumPart = actualNumPart;
    return simu;
  }
};

SimulationLoader *multidarkLoader(const string& multidarkname)
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

  return new MultiDarkLoader(multidarkname, header);
}

  

