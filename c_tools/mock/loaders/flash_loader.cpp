/*+
    VIDE -- Void IDentification and Examination -- ./c_tools/mock/loaders/flash_loader.cpp
    Copyright (C) 2010-2013 Guilhem Lavaux
    Copyright (C) 2011-2013 P. M. Sutter

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
#include <CosmoTool/loadFlash.hpp>
#include <CosmoTool/fortran.hpp>
#include "simulation_loader.hpp"

using namespace std;
using namespace CosmoTool;

class FlashLoader: public SimulationLoader
{
private:
  int load_flags;
  bool onefile;
  int _num_files;
  SimuData *gadget_header;
  string snapshot_name;
  SimulationPreprocessor *preproc;
public:
  FlashLoader(const string& basename, SimuData *header, int flags, bool singleFile, int _num, SimulationPreprocessor *p)
    : snapshot_name(basename), load_flags(flags), onefile(singleFile), _num_files(_num), gadget_header(header), preproc(p)
  {
  }
  
  ~FlashLoader()
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
      d = loadFlashMulti(snapshot_name.c_str(), -1, load_flags);
    else
      d = loadFlashMulti(snapshot_name.c_str(), id, load_flags);

    if (d->Id != 0)
      {
	long *uniqueID = new long[d->NumPart];
	for (long i = 0; i < d->NumPart; i++)
	  {
	    uniqueID[i] = d->Id[i];
	  }
	d->new_attribute("uniqueID", uniqueID, delete_adaptor<long>);
      }

    applyTransformations(d);
    basicPreprocessing(d, preproc);

    return d;
  }
};


SimulationLoader *flashLoader(const std::string& snapshot, int flags, SimulationPreprocessor *p)
{
  bool singleFile;
  int num_files;
  SimuData *d;

  try
    {
      d = loadFlashMulti(snapshot.c_str(), -1, 0);
      singleFile = true;
      num_files = 1;
    }
  catch (const NoSuchFileException& e)
    {
      try
        {
          d = loadFlashMulti(snapshot.c_str(), 0, 0);
	  num_files = 0;
        }
      catch(const NoSuchFileException& e)
        {
          return 0;
        }
    }
    
  assert(d != 0);
  SimuData *header = d;
    
  if (!singleFile)
    {
      try
	{
	  while ((d = loadFlashMulti(snapshot.c_str(), num_files, 0)) != 0)
	    {
	      num_files++;
	      delete d;
	    }
	}
      catch(const NoSuchFileException& e)
	{
	}
    }
    
  return new FlashLoader(snapshot, header, flags, singleFile, num_files, p);
}
