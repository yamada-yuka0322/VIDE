/*+
    VIDE -- Void IDentification and Examination -- ./c_tools/mock/loaders/gadget_loader.cpp
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



#include <vector>
#include <cassert>
#include <string>
#include <CosmoTool/loadGadget.hpp>
#include <CosmoTool/fortran.hpp>
#include "simulation_loader.hpp"

using namespace std;
using namespace CosmoTool;

class GadgetLoader: public SimulationLoader
{
private:
  int load_flags;
  bool onefile;
  int _num_files;
  double unitMpc;
  int gadgetFormat;
  SimuData *gadget_header;
  string snapshot_name;
  SimulationPreprocessor *preproc;
public:
  GadgetLoader(const string& basename, SimuData *header, int flags, bool singleFile, int _num, double unit, int gadgetFormat, SimulationPreprocessor *p)
    : snapshot_name(basename), load_flags(flags), onefile(singleFile), _num_files(_num), unitMpc(1/unit), gadget_header(header), gadgetFormat(gadgetFormat), preproc(p)
  {
  }
  
  ~GadgetLoader()
  {
    delete gadget_header;
  }
  
  SimuData *getHeader() {
    return gadget_header;
  }
  
  int num_files() {
    return _num_files;
  }

  SimuData *loadFile(int id) {
    SimuData *d;
    
    if (onefile && id > 0)
      return 0;
    if (id >= _num_files)
      return 0;

    if (onefile)
      d = loadGadgetMulti(snapshot_name.c_str(), -1, load_flags, gadgetFormat);
    else
      d = loadGadgetMulti(snapshot_name.c_str(), id, load_flags, gadgetFormat);

    if (d->Id != 0)
      {
	long *uniqueID = new long[d->NumPart];
	for (long i = 0; i < d->NumPart; i++)
	  {
	    uniqueID[i] = d->Id[i];
	  }
	d->new_attribute("uniqueID", uniqueID, delete_adaptor<long>);
      }

    for (int k = 0; k < 3; k++)
      {
        if (d->Pos[k] != 0)
          {
            for (long i = 0; i < d->NumPart; i++)
              d->Pos[k][i] *= unitMpc;
          }
      }
    d->BoxSize *= unitMpc;

    applyTransformations(d);
    basicPreprocessing(d, preproc);

    return d;
  }
};


SimulationLoader *gadgetLoader(const std::string& snapshot, double Mpc_unitLength, int flags, int gadgetFormat, SimulationPreprocessor *p)
{
  bool singleFile = false;
  int num_files;
  SimuData *d;

  cout << " File to load is:" << snapshot.c_str() << endl;

  try
    {
      d = loadGadgetMulti(snapshot.c_str(), -1, 0, gadgetFormat);
      singleFile = true;
      num_files = 1;
    }
  catch (const NoSuchFileException& e)
    {
      try
        {
          d = loadGadgetMulti(snapshot.c_str(), 0, 0, gadgetFormat);
	  num_files = 0;
        }
      catch(const NoSuchFileException& e)
        {
          return 0;
        }
    }
    
  assert(d != 0);
  SimuData *header = d;

  header->BoxSize /= Mpc_unitLength;
    
  if (!singleFile)
    {
      try
	{
	  while ((d = loadGadgetMulti(snapshot.c_str(), num_files, 0, gadgetFormat)) != 0)
	    {
	      num_files++;
	      delete d;
	    }
	}
      catch(const NoSuchFileException& e)
	{
	}
    }
    
  return new GadgetLoader(snapshot, header, flags, singleFile, num_files, Mpc_unitLength, gadgetFormat, p);
}
