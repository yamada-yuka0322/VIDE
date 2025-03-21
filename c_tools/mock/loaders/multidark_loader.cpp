/*+
    VIDE -- Void IDentification and Examination -- ./c_tools/mock/loaders/multidark_loader.cpp
    Copyright (C) 2010-2014 Guilhem Lavaux
    Copyright (C) 2011-2014 P. M. Sutter

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; version 2 of the License.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
+*/



#include <cassert>
#include <boost/format.hpp>
#include <string>
#include <fstream>
#include <iostream>
#include <CosmoTool/loadRamses.hpp>
#include <CosmoTool/fortran.hpp>
#include "simulation_loader.hpp"

using boost::format;
using namespace std;
using namespace CosmoTool;

static const double one_kpc = 3.08567802e16; /* km */
static const double one_Gyr = 3.1558149984e16; /* sec */

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

    for (int k = 0; k < 3; k++) {
      simu->Pos[k] = new float[allocated];
      simu->Vel[k] = new float[allocated];
    }
    simu->Id = new int64_t[allocated];
    long *uniqueID = new long[allocated];
    long *index = new long[allocated];

    double tempData;

    double shift = 0.5*simu->BoxSize;
    double rescale_position = simu->Hubble*1.e-5*100./simu->time;
    double rescale_velocity = one_kpc/one_Gyr;
#define INFINITY std::numeric_limits<float>::max()
    float min_pos[3] = {INFINITY,INFINITY, INFINITY}, max_pos[3] = {-INFINITY,-INFINITY,-INFINITY};

    cout << "loading multidark particles" << endl;
    long actualNumPart = 0;

    for (long i = 0; i < simu->NumPart; i++) {
      SingleParticle p;

      fp >> p.ID >> p.Pos[0] >> p.Pos[1]
         >> p.Pos[2] >> p.Vel[0] >> p.Vel[1] >> p.Vel[2] >> tempData;

      if (p.ID == -99 && 
          p.Pos[0] == -99 && p.Pos[1] == -99 && 
          p.Pos[2] == -99 && p.Vel[2] == -99)
        break;

      //p.Pos[0] = p.Pos[0]*rescale_position + shift;
      //p.Pos[1] = p.Pos[1]*rescale_position + shift;
      //p.Pos[2] = p.Pos[2]*rescale_position + shift;
      //p.Vel[2] = p.Vel[2]*rescale_velocity;

      // enforce box size in case of roundoff error
      for (int k = 0; k < 3; k++) {
        if (p.Pos[k] < 0) p.Pos[k] += simu->BoxSize;
        if (p.Pos[k] >= simu->BoxSize) p.Pos[k] -= simu->BoxSize;
      }

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

       for (int k = 0; k < 3; k++) {
         min_pos[k] = std::min(min_pos[k], p.Pos[k]);
         max_pos[k] = std::max(max_pos[k], p.Pos[k]);
       }
    }
    for (int k = 0; k < 3; k++) cout << boost::format("min[%d] = %g, max[%d] = %g") % k % min_pos[k] % k %max_pos[k] << endl;

    applyTransformations(simu);
    simu->NumPart = actualNumPart;
    simu->TotalNumPart = actualNumPart;
    simu->new_attribute("uniqueID", uniqueID, delete_adaptor<long>);
    simu->new_attribute("index", index, delete_adaptor<long>);
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

  

