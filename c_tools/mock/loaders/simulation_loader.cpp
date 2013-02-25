#include <cmath>
#include <CosmoTool/loadSimu.hpp>
#include "simulation_loader.hpp"

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
