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
public:
  RamsesLoader(const string& basename, int baseid, bool dp, SimuData *header, int flags, int _num)
    : snapshot_name(basename), load_flags(flags), _num_files(_num), double_precision(dp),
      ramses_header(header)
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
      
    long *uniqueID = new long[d->NumPart];
    for (long i = 0; i < d->NumPart; i++)
      {
        uniqueID[i] = d->Id[i];
      }

    d->new_attribute("uniqueID", uniqueID, delete_adaptor<long>);

    applyTransformations(d);

    return d;
  }
};

SimulationLoader *ramsesLoader(const std::string& snapshot, int baseid, bool double_precision, int flags)
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

  return new RamsesLoader(snapshot, baseid, double_precision, header, flags, num_files);
}

