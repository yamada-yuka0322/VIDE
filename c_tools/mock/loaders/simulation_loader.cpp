#include <cmath>
#include <CosmoTool/loadSimu.hpp>
#include "simulation_loader.hpp"
#ifdef SDF_SUPPORT
#include "sdfloader_internal.hpp"
#endif

using std::min;
using namespace CosmoTool;

template<typename T> void reallocArray(T *& a, long newSize, long toCopy)
{
  T *b = new T[newSize];
  if (a != 0)
    {
      memcpy(b, a, sizeof(T)*toCopy);
      delete[] a;
    }
  a = b;
}

void SimulationLoader::reallocSimu(SimuData *s, long newNumPart)
{
  long to_copy = min(newNumPart, s->NumPart);
  
  for (int j = 0; j < 3; j++)
    {
      reallocArray(s->Pos[j], newNumPart, to_copy);
      reallocArray(s->Vel[j], newNumPart, to_copy);
    }
  reallocArray(s->Id, newNumPart, to_copy);
  reallocArray(s->type, newNumPart, to_copy);
}

void SimulationLoader::applyTransformations(SimuData *s)
{
  float redshift_gravity = do_redshift ? 1.0 : 0.0;

  for (int i = 0; i < s->NumPart; i++)
    {
      s->Pos[redshift_axis][i] += 
	redshift_gravity*s->Vel[redshift_axis][i]/100.;
    }
}


void SimulationLoader::basicPreprocessing(SimuData *d, 
                                          SimulationPreprocessor *preproc)
{
  if (preproc == 0)
    return;

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
}

void simulationLoadersInit(int& argc, char **& argv)
{
#ifdef SDF_SUPPORT
  sdfLoaderInit(argc, argv);
#endif
}
