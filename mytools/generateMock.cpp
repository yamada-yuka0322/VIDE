#include <cassert>
#include <iostream>
#include <fstream>
#include <CosmoTool/loadSimu.hpp>
#include <CosmoTool/loadRamses.hpp>
#include "generateMock_conf.h"

using namespace std;
using namespace CosmoTool;

SimuData *doLoadRamses(const char *basename, int baseid, int velAxis, bool goRedshift)
{
  SimuData *d, *outd;

  d = loadRamsesSimu(basename, baseid, -1, 0);
  outd = new SimuData;

  outd->NumPart = d->TotalNumPart;
  outd->BoxSize = d->BoxSize;
  outd->TotalNumPart = outd->NumPart;
  outd->Hubble = d->Hubble;
  for (int k = 0; k < 3; k++)
    outd->Pos[k] = new float[outd->NumPart];
  outd->Vel[2] = new float[outd->NumPart];
  delete d;

  int curCpu = 0;
  cout << "loading cpu 0 " << endl;
  while (d = loadRamsesSimu(basename, baseid, curCpu, NEED_POSITION|NEED_VELOCITY|NEED_GADGET_ID))
    { 
       for (int k = 0; k < 3; k++)
         for (int i = 0; i < d->NumPart; i++)
           {
             assert(d->Id[i] >= 1);
             assert(d->Id[i] <= outd->TotalNumPart);
             outd->Pos[k][d->Id[i]-1] = d->Pos[k][i];
             outd->Vel[2][d->Id[i]-1] = d->Vel[velAxis][i];
           }

       if (goRedshift)
         for (int i = 0; i < d->NumPart; i++)
            outd->Pos[velAxis][d->Id[i]-1] += d->Vel[velAxis][i]/100.;

       delete d;
       curCpu++;
       cout << "loading cpu " << curCpu  << endl;
    }

  return outd;
}

int main(int argc, char **argv)
{
  generateMock_info args_info;
  generateMock_conf_params args_params;
  SimuData *simu;
 
  generateMock_conf_init(&args_info);
  generateMock_conf_params_init(&args_params);
  
  args_params.check_required = 0;
  if (generateMock_conf_ext (argc, argv, &args_info, &args_params))
    return 1;
  
  if (!args_info.configFile_given)
    {
      if (generateMock_conf_required (&args_info, GENERATEMOCK_CONF_PACKAGE))
        return 1;
    }
  else
    {
      args_params.check_required = 1;
      args_params.initialize = 0;
      if (generateMock_conf_config_file (args_info.configFile_arg,
					 &args_info,
					 &args_params))
	return 1;
    }
  
  generateMock_conf_print_version();
  
  simu = doLoadRamses(args_info.ramsesBase_arg, 
		      args_info.ramsesId_arg, 
		      args_info.axis_arg, true);

  cout << "Hubble = " << simu->Hubble << endl;
  cout << "Boxsize = " << simu->BoxSize << endl;
  
  return 0;
}
