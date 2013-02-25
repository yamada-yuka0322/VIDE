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
public:
  FlashLoader(const string& basename, SimuData *header, int flags, bool singleFile, int _num)
    : snapshot_name(basename), load_flags(flags), onefile(singleFile), _num_files(_num), gadget_header(header)
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
    
  return new FlashLoader(snapshot, header, flags, singleFile, num_files);
}
