/*+
    VIDE -- Void IDEntification pipeline -- ./c_tools/mock/loaders/basic_loader.cpp
    Copyright (C) 2010-2013 Guilhem Lavaux
    Copyright (C) 2011-2013 Paul M. Sutter

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

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
#include <iostream>
#include <boost/format.hpp>
#include <fstream>
#include <CosmoTool/yorick.hpp>
#include "simulation_loader.hpp"

using namespace std;
using namespace CosmoTool;
using boost::format;
using boost::str;

class BasicGroupLoader: public SimulationLoader
{
private:
  string storage;
  SimuData *header;
  int flags, numFiles;
public:
  BasicGroupLoader(const string& storage_path, SimuData *header, int flags, int numfiles)
  {
    this->header = header;
    this->storage = storage_path;
    this->flags = flags;
    this->numFiles = numfiles;
  }
  
  SimuData *getHeader() {
    return header;
  }
  
  int num_files() {
    return numFiles;
  }
  
  SimuData *loadFile(int id) {    
    if (id < 0 || id >= numFiles)
      return 0;

    SimuData *simu = new SimuData;
    uint32_t *dimlist, dimrank;
    string fname;

    simu->time = header->time;
    simu->TotalNumPart = header->NumPart;
    simu->Omega_M = header->Omega_M;
    simu->Omega_Lambda = header->Omega_Lambda;
    simu->NumPart = -1;

    if (flags & NEED_POSITION)
      {
	loadArray(str(format("%s/x_%d.nc") % storage % id), 
		  simu->Pos[0], dimlist, dimrank);
	assert(dimrank == 1);
	simu->NumPart = dimlist[0];
	
	loadArray(str(format("%s/y_%d.nc") % storage % id), 
		  simu->Pos[1], dimlist, dimrank);
	assert(dimrank == 1);
	assert(simu->NumPart == dimlist[0]);
	
	loadArray(str(format("%s/z_%d.nc") % storage % id), 
		  simu->Pos[2], dimlist, dimrank);
	assert(dimrank == 1);
	assert(simu->NumPart == dimlist[0]);
      }

    if (flags & NEED_VELOCITY)
      {
	loadArray(str(format("%s/vx_%d.nc") % storage % id), 
		  simu->Vel[0], dimlist, dimrank);
	assert(dimrank == 1);
	if (simu->NumPart < 0)
	  simu->NumPart = dimlist[0];

	assert(simu->NumPart == dimlist[0]);
	
	loadArray(str(format("%s/vy_%d.nc") % storage % id), 
		  simu->Vel[0], dimlist, dimrank);
	assert(dimrank == 1);
	assert(simu->NumPart == dimlist[0]);
	
	loadArray(str(format("%s/vz_%d.nc") % storage % id), 
		  simu->Vel[2], dimlist, dimrank);
	assert(dimrank == 1);
	assert(simu->NumPart == dimlist[0]);
      }

    if (flags & NEED_GADGET_ID)
      {
	loadArray(str(format("%s/id_%d.nc") % storage % id), 
		  simu->Id, dimlist, dimrank);
	assert(dimrank == 1);
	if (simu->NumPart < 0)
	  simu->NumPart = dimlist[0];

	assert(simu->NumPart == dimlist[0]);
      }
    
    return simu;
  }

  ~BasicGroupLoader()
  {
    delete header;
  }
};

SimulationLoader *basicGroupLoader(const std::string& simupath, 
				   int flags)
{
  SimuData *header;
  ifstream f;
  string header_path = simupath + "/header.txt";
  int numFiles;

  header = new SimuData;
  f.open(header_path.c_str());

  f >> header->time
    >> header->Omega_M
    >> header->Omega_Lambda 
    >> header->NumPart
    >> numFiles;

  return new BasicGroupLoader(simupath, header, flags, numFiles);
}
