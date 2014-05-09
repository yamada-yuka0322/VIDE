/*+
    VIDE -- Void IDentification and Examination -- ./c_tools/libzobov/loadZobov.cpp
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
#include <string>
#include <fstream>
#include <iostream>
#include <vector>
#include <cstdlib>
#include <sstream>
#include <algorithm>
#include "loadZobov.hpp"
#include <CosmoTool/fortran.hpp>

using namespace std;
using namespace CosmoTool;

bool loadZobov(const char *descName, const char *adjName, const char *voidsName, 
	       const char *volName, ZobovRep& z)
{
  ifstream descFile(descName);
  ifstream adjFile(adjName);
  ifstream volFile(voidsName);
  int32_t numParticles, numZones, numPinZone;
  int32_t totalParticles;
  int32_t numVoids;
  int32_t minParticlesInZone, maxParticlesInZone;

  adjFile.read((char *)&numParticles, sizeof(numParticles));
  adjFile.read((char *)&numZones, sizeof(numZones));
  if (!adjFile)
    return false;

  cout << "Number of particles = " << numParticles << endl;
  cout << "Number of zones = " << numZones << endl;

  totalParticles = 0;

  minParticlesInZone = -1;
  maxParticlesInZone = -1;

  z.allZones.resize(numZones);
  for (int zone = 0; zone < numZones; zone++)
    {
      adjFile.read((char *)&numPinZone, sizeof(numPinZone));
      if (!adjFile)
	{
	  cout << "Problem on the zone " << zone << " / " << numZones << endl;
	  return false;
	}
      
      z.allZones[zone].pId.resize(numPinZone);
      adjFile.read((char *)&z.allZones[zone].pId[0], sizeof(int)*numPinZone);

      if (maxParticlesInZone < 0 || numPinZone > maxParticlesInZone)
	maxParticlesInZone = numPinZone;

      if (minParticlesInZone < 0 || numPinZone < minParticlesInZone)
	minParticlesInZone = numPinZone;

      totalParticles += numPinZone;      
    }
  cout << "Zoned " << totalParticles << endl;

  cout << "Minimum number of particles in zone = " << minParticlesInZone << endl;
  cout << "Maximum number of particles in zone = " << maxParticlesInZone << endl;
  
  if (totalParticles != numParticles)
    {
      cerr << "The numbers of particles are inconsistent ! (" << totalParticles << " vs " << numParticles << ")"<< endl;
      abort();
    }

  volFile.read((char *)&numVoids, sizeof(numVoids));
  if (!volFile)
    return false;

  cout << "Number of voids = " << numVoids << endl;

  z.allVoids.resize(numVoids);
  for (int v = 0; v < numVoids; v++)
    {
      int32_t numZinV;

      volFile.read((char *)&numZinV, sizeof(numZinV));
      if (!volFile)
	return false;

      z.allVoids[v].zId.resize(numZinV);

      int *zId = new int[numZinV];      

      volFile.read((char *)zId, sizeof(int)*numZinV);
      for (int k = 0; k < numZinV; k++)
	z.allVoids[v].zId[k] = zId[k];
      std::sort(&z.allVoids[v].zId[0], &z.allVoids[v].zId[numZinV]);

      delete[] zId;
    }

  if (volName != 0)
    {
      cout << "Loading particle volumes (requested)" << endl;
      ifstream f(volName);
      int numParticles;

      if (!f)
	{
	  cerr << "No such file " << volName << endl;
	  abort();
	}

      f.read((char *)&numParticles, sizeof(int));
      z.particleVolume.resize(numParticles);
      f.read((char *)&z.particleVolume[0], sizeof(float)*numParticles);
    }

  cout << "Loading description" << endl;
  int lineid = 1;

  string line;
  getline(descFile, line);
  lineid++;
  getline(descFile, line);
  lineid++;
  getline(descFile, line);
  lineid++;
  while (!descFile.eof())
    {
      istringstream lineStream(line.c_str());
      int orderId, volId, coreParticle, numParticlesInZone, numZonesInVoid, numInVoid;
      float coreDensity, volumeZone, volumeVoid, densityContrast;
      float probability;

      lineStream
	>> orderId
	   >> volId
	   >> coreParticle
	   >> coreDensity
	   >> volumeZone
	   >> numParticlesInZone
	   >> numZonesInVoid
	   >> volumeVoid
	   >> numInVoid
	   >> densityContrast
	   >> probability;
      if (!lineStream)
	{
	  cerr << "Error in text stream at line " << lineid << endl;
	  abort();
	}
    
/*  
    sscanf(line.c_str(), "%d %d %d %f %f %d %d %f %d %f %f\n", &orderId, &volId,
           &coreParticle, &coreDensity, &volumeZone, &numParticlesInZone, 
           &numZonesInVoid,
           &volumeVoid, &numInVoid, &densityContrast, &probability);
*/

      z.allVoids[volId].proba = probability;
      z.allVoids[volId].volume = volumeVoid;
      z.allVoids[volId].numParticles = numInVoid;
      z.allVoids[volId].coreParticle = coreParticle;
      
      // Sanity check
      int actualNumber = 0;
      for (int j = 0; j < z.allVoids[volId].zId.size(); j++)
        {
          int zzid = z.allVoids[volId].zId[j];
          assert(zzid < z.allZones.size());
	  actualNumber += z.allZones[zzid].pId.size();
        }
      
      if (actualNumber != numInVoid)
	{
	  cerr << "Sanity check failed."
	       << " The number of particles in the description (" 
	       << numInVoid
	       << ") is different from the one in the file ("
	       << actualNumber << ")" << endl;
	}
      getline(descFile, line);
      lineid++;
    }

  cout << "Done loading" << endl;


  

  return true;
}

bool loadZobovParticles(const char *fname, std::vector<ZobovParticle>& particles)
{
  UnformattedRead f(fname);
  int N;

  f.beginCheckpoint();
  N = f.readInt32();
  f.endCheckpoint();

  particles.resize(N);

  f.beginCheckpoint();
  for (int i = 0; i < N; i++)
    {
      particles[i].x = f.readReal32();
    }
  f.endCheckpoint();

  f.beginCheckpoint();
  for (int i = 0; i < N; i++)
    {
      particles[i].y = f.readReal32();
    }
  f.endCheckpoint();

  f.beginCheckpoint();
  for (int i = 0; i < N; i++)
    {
      particles[i].z = f.readReal32();
    }
  f.endCheckpoint();

  return true;
}
