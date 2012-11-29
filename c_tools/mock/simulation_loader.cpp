#include <CosmoTool/loadSimu.hpp>
#include "simulation_loader.hpp"

using namespace CosmoTool;

void SimulationLoader::applyTransformations(SimuData *s)
{
  float redshift_gravity = do_redshift ? 1.0 : 0.0;

  for (int i = 0; i < s->NumPart; i++)
    {
      s->Pos[redshift_axis][i] += 
	redshift_gravity*s->Vel[redshift_axis][i]/100.;
    }
}
