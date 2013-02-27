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
  SimuData *gadget_header;
  string snapshot_name;
  SimulationPreprocessor *preproc;
public:
  GadgetLoader(const string& basename, SimuData *header, int flags, bool singleFile, int _num, double unit, SimulationPreprocessor *p)
    : snapshot_name(basename), load_flags(flags), onefile(singleFile), _num_files(_num), unitMpc(1/unit), gadget_header(header), preproc(p)
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
      d = loadGadgetMulti(snapshot_name.c_str(), -1, load_flags);
    else
      d = loadGadgetMulti(snapshot_name.c_str(), id, load_flags);

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

    long numAccepted = 0;
    bool *accepted = new bool[d->NumPart];
    for (long i = 0; i < d->NumPart; i++)
      {
	SingleParticle p;

	for (int k = 0; k < 3; k++)
	  {
	    p.Pos[k] = (d->Pos[k]) ? 0 : d->Pos[k][i];
	    p.Vel[k] = (d->Vel[k]) ? 0 : d->Vel[k][i];
	  }
	p.ID = (d->Id) ? 0 : d->Id[i];
	
	accepted[i] = preproc->accept(p);
	numAccepted += accepted[i];
      }
    
    for (int k = 0; k < 3; k++)
      {
	filteredCopy(d->Pos[k], accepted, d->NumPart);
	filteredCopy(d->Vel[k], accepted, d->NumPart);
      }
    filteredCopy(d->Id, accepted, d->NumPart);
    filteredCopy(d->type, accepted, d->NumPart);
    delete[] accepted;
    d->NumPart = numAccepted;

    return d;
  }
};


SimulationLoader *gadgetLoader(const std::string& snapshot, double Mpc_unitLength, int flags, SimulationPreprocessor *p)
{
  bool singleFile = false;
  int num_files;
  SimuData *d;

  try
    {
      d = loadGadgetMulti(snapshot.c_str(), -1, 0);
      singleFile = true;
      num_files = 1;
    }
  catch (const NoSuchFileException& e)
    {
      try
        {
          d = loadGadgetMulti(snapshot.c_str(), 0, 0);
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
	  while ((d = loadGadgetMulti(snapshot.c_str(), num_files, 0)) != 0)
	    {
	      num_files++;
	      delete d;
	    }
	}
      catch(const NoSuchFileException& e)
	{
	}
    }
    
  return new GadgetLoader(snapshot, header, flags, singleFile, num_files, Mpc_unitLength, p);
}
