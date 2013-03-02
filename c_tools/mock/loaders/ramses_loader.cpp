/*+
    VIDE -- Void IDEntification pipeline -- ./c_tools/mock/loaders/ramses_loader.cpp
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
#include <CosmoTool/loadRamses.hpp>
#include <CosmoTool/fortran.hpp>
#include "simulation_loader.hpp"

using namespace std;
using namespace CosmoTool;

class RamsesLoader: public SimulationLoader
{
private:
  int load_flags;
  int _num_files;
  int baseid;
  bool double_precision;
  SimuData *ramses_header;
  string snapshot_name;
  SimulationPreprocessor *preproc;
public:
  RamsesLoader(const string& basename, int baseid, bool dp, SimuData *header, int flags, int _num, SimulationPreprocessor *p)
    : snapshot_name(basename), load_flags(flags), _num_files(_num), double_precision(dp),
      ramses_header(header), preproc(p)
  {
  }
  
  ~RamsesLoader()
  {
    delete ramses_header;
  }
  
  SimuData *getHeader() {
    return ramses_header;
  }
  
  int num_files() {
    return _num_files;
  }
  
  SimuData *loadFile(int id) {
    SimuData *d;
    
    if (id >= _num_files)
      return 0;

    d = loadRamsesSimu(snapshot_name.c_str(), baseid, id, double_precision, load_flags);
    assert(d != 0);

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

SimulationLoader *ramsesLoader(const std::string& snapshot, int baseid, bool double_precision, int flags, SimulationPreprocessor *p)
{
  SimuData *d, *header;
  int num_files = 0;

  header = loadRamsesSimu(snapshot.c_str(), baseid, 0, double_precision, 0);
  if (header == 0)
    return 0;

  while ((d = loadRamsesSimu(snapshot.c_str(), baseid, num_files, double_precision, 0)) != 0)
    {
      num_files++;
      delete d;
    }

  return new RamsesLoader(snapshot, baseid, double_precision, header, flags, num_files, p);
}

