/*+
    VIDE -- Void IDentification and Examination -- ./c_tools/mock/loaders/ramses_loader.cpp
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
#include <CosmoTool/loadRamses.hpp>
#include <CosmoTool/fortran.hpp>
#include "simulation_loader.hpp"

// ben edit

#include <iostream>
#include <iomanip>
#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <sys/types.h>
#include <regex.h>
#include <cstring>
#include <cstdlib>
#include <cassert>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <map>


// end ben edit


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
    

// cludgy hack until I can get baseid working in this function... the problem is with   SimuData *loadFile(int id) ... just need to learn a bit more C++ ~ Ben

    string baseidstr = snapshot_name.c_str();
    unsigned found = baseidstr.find_last_of("/");
    baseidstr = baseidstr.substr(found-5,found);
    baseidstr = baseidstr.substr(0,5);  // remove trailing slash
    baseid = atoi(baseidstr.c_str());
    
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
  int num_files = 0; // how many particle files are there?

  header = loadRamsesSimu(snapshot.c_str(), baseid, 0, double_precision, 0);

  if (header == 0)
    return 0;


// count number of CPU's for this output... kinda copy/paste of loadRamses.cpp in cosmotool/src

  ostringstream ss_fname;  
  ss_fname << snapshot.c_str() << "/info_" << setfill('0') << setw(5) << baseid << ".txt";
  cout << "Opening info file " << ss_fname.str() << " to find cpu number" << endl;
  ifstream infile(ss_fname.str().c_str());
  if (!infile)
    return 0;


  int err;
  regex_t unit_l_rx;

  //  const char *pattern = "^unit_l[ ]*=[ ]*([0-9\\.E+\\-]+)";
  const char *pattern = "^([A-Za-z_]+)[ ]*=[ ]*([0-9\\.E+\\-]+)";

  err = regcomp (&unit_l_rx, pattern, REG_EXTENDED);
  cout << unit_l_rx.re_nsub << endl;
  if (err)
    {
      char errString[255];
      regerror(err, &unit_l_rx, errString, sizeof(errString));
      cout << errString << endl;
      abort();
    }

  map<string,double> infoMap;
  string line;
  while (getline(infile, line))
    {
      regmatch_t allMatch[4];
      if (!regexec(&unit_l_rx, line.c_str(), 4, allMatch, 0))
	{
	  uint32_t start0 = allMatch[1].rm_so, end0 = allMatch[1].rm_eo;
	  uint32_t start1 = allMatch[2].rm_so, end1 = allMatch[2].rm_eo;

	  string keyword = line.substr(start0, end0-start0);
	  istringstream iss(line.substr(start1, end1-start1));
	  double unitLength;

	  iss >> unitLength;

	  infoMap[keyword] = unitLength;
	}
    }
  
  regfree(&unit_l_rx);


  num_files = infoMap["ncpu"];


  return new RamsesLoader(snapshot, baseid, double_precision, header, flags, num_files, p);
}

