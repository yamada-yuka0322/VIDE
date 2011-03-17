#include <cmath>
#include <cassert>
#include <iostream>
#include <fstream>
#include <string>
#include <CosmoTool/loadSimu.hpp>
#include <CosmoTool/loadRamses.hpp>
#include <CosmoTool/loadGadget.hpp>
#include <CosmoTool/interpolate.hpp>
#include <CosmoTool/fortran.hpp>
#include "generateMock_conf.h"
#include "gslIntegrate.hpp"
#include <netcdfcpp.h>

using namespace std;
using namespace CosmoTool;

#define LIGHT_SPEED 299792.458

SimuData *doLoadRamses(const char *basename, int baseid, int velAxis, bool goRedshift)
{
  SimuData *d, *outd;

  d = loadRamsesSimu(basename, baseid, -1, 0);
  outd = new SimuData;

  outd->NumPart = d->TotalNumPart;
  outd->BoxSize = d->BoxSize;
  outd->TotalNumPart = outd->NumPart;
  outd->Hubble = d->Hubble;
  outd->Omega_Lambda = d->Omega_Lambda;
  outd->Omega_M = d->Omega_M;
  outd->time = d->time;

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

SimuData *doLoadGadget(const char *gadgetname, int velAxis, bool goRedshift)
{
  SimuData *d, *outd;

  try
    {
      d = loadGadgetMulti(gadgetname, -1, 0);
    }
  catch (const NoSuchFileException& e)
    {
      try
	{
	  d = loadGadgetMulti(gadgetname, 0, 0);
	}
      catch(const NoSuchFileException& e)
	{
	  return 0;
	}
    }
  outd = new SimuData;

  outd->NumPart = d->TotalNumPart;
  outd->BoxSize = d->BoxSize/1000;
  outd->TotalNumPart = outd->NumPart;
  outd->Hubble = d->Hubble;
  outd->Omega_Lambda = d->Omega_Lambda;
  outd->Omega_M = d->Omega_M;
  outd->time = d->time;

  for (int k = 0; k < 3; k++)
    outd->Pos[k] = new float[outd->NumPart];
  outd->Vel[2] = new float[outd->NumPart];
  delete d;

  int curCpu = 0;
  cout << "loading file 0 " << endl;
  try
    {
      while (1)
	{
	  d = loadGadgetMulti(gadgetname, curCpu, NEED_POSITION|NEED_VELOCITY|NEED_GADGET_ID);
	  for (int k = 0; k < 3; k++)
	    for (int i = 0; i < d->NumPart; i++)
	      {
		assert(d->Id[i] >= 1);
		assert(d->Id[i] <= outd->TotalNumPart);
		outd->Pos[k][d->Id[i]-1] = d->Pos[k][i]/1000;
		outd->Vel[2][d->Id[i]-1] = d->Vel[velAxis][i];
	      }
	  
	  if (goRedshift)
	    for (int i = 0; i < d->NumPart; i++)
	      outd->Pos[velAxis][d->Id[i]-1] += d->Vel[velAxis][i]/100.;

	  delete d;
	  curCpu++;
	  cout << "loading file " << curCpu  << endl;
	}
    }
  catch (const NoSuchFileException& e)
    {
    }

  return outd;
}

static double cubic(double a)
{
  return a*a*a;
}

struct TotalExpansion
{
  double Omega_M, Omega_L;

  double operator()(double z)
  {
    return 1/sqrt(Omega_M*cubic(1+z) + Omega_L);
  }
};

Interpolate make_cosmological_redshift(double OM, double OL, double z0, double z1, int N = 1000)
{
  TotalExpansion e_computer;
  double D_tilde, Q, Qprime;
  InterpolatePairs pairs;
  
  e_computer.Omega_M = OM;
  e_computer.Omega_L = OL;
  
  pairs.resize(N);
  ofstream f("comoving_distance.txt");

  for (int i = 0; i < N; i++)
    {      
      double z = z0 + (z1-z0)/N*i;

      pairs[i].second = z;
      pairs[i].first = gslIntegrate(e_computer, 0, z, 1e-3);
      f << z << " " <<  pairs[i].first << endl;
    }

  return buildFromVector(pairs);
}

void metricTransform(SimuData *data, int axis, bool reshift, bool pecvel, double*& expfact)
{
  int x0, x1, x2;

  switch (axis) {
  case 0:
    x0 = 1; x1 = 2; x2 = 0;
    break;
  case 1:
    x0 = 0; x1 = 2; x2 = 1;
    break;
  case 2:
    x0 = 0; x1 = 1; x2 = 2;
    break;
  default:
    abort();
  }

  Interpolate z_vs_D = make_cosmological_redshift(data->Omega_M, data->Omega_Lambda, 0., 2.0); // Redshift 2 should be sufficient ?
    
  double z0 = 1/data->time - 1;
  double z_base = reshift ? z0 : 0;
  TotalExpansion e_computer;
  double baseComovingDistance;

  expfact = new double[data->NumPart];

  cout << "Using base redshift z=" << z0 << endl;

  e_computer.Omega_M = data->Omega_M;
  e_computer.Omega_L = data->Omega_Lambda;
  baseComovingDistance = LIGHT_SPEED/100.* gslIntegrate(e_computer, 0, z0, 1e-3);
  cout << "Comoving distance = " << baseComovingDistance << " Mpc/h" << endl;

  for (uint32_t i = 0; i < data->NumPart; i++)
    {
      float& x = data->Pos[x0][i];
      float& y = data->Pos[x1][i];
      float& z = data->Pos[x2][i];
      float& v = data->Vel[2][i];
      float z_old = z;
      
      double reduced_red = (z + baseComovingDistance)*100./LIGHT_SPEED;
      try
        {

          // Distorted redshift
          if (reduced_red == 0)
            z = 0;
          else
            z = (z_vs_D.compute(reduced_red)-z_base)*LIGHT_SPEED/100.;
          expfact[i] = z / z_old; 
         // Add peculiar velocity
         if (pecvel)
	  z += v/100;
       }
      catch(const InvalidRangeException& e) {
       cout << "Trying to interpolate out of the tabulated range." << endl;
       cout << "The offending value is z=" << reduced_red << endl;
       abort();
      }
    }
}

void generateOutput(SimuData *data, int axis, 
		    const std::string& fname)
{
  UnformattedWrite f(fname);

  cout << "Generating output particles to " << fname << endl;
 
  int x0, x1, x2;
  
  switch (axis) {
  case 0:
    x0 = 1; x1 = 2; x2 = 0;
    break;
  case 1:
    x0 = 0; x1 = 2; x2 = 1;
    break;
  case 2:
    x0 = 0; x1 = 1; x2 = 2;
    break;
  default:
    abort();
  }

  f.beginCheckpoint();
  f.writeInt32(data->NumPart);
  f.endCheckpoint();

  cout << "Writing X components..." << endl;
  f.beginCheckpoint();
  for (uint32_t i = 0; i < data->NumPart; i++)
    {
      f.writeReal32(data->Pos[x0][i]);
    }
  f.endCheckpoint();

  cout << "Writing Y components..." << endl;
  f.beginCheckpoint();
  for (uint32_t i = 0; i < data->NumPart; i++)
    {
      f.writeReal32(data->Pos[x1][i]);
    }
  f.endCheckpoint();
  
  cout << "Writing Z components..." << endl;
  f.beginCheckpoint();
  for (uint32_t i = 0; i < data->NumPart; i++)
    {
      f.writeReal32(data->Pos[x2][i]);
    }
  f.endCheckpoint();
}

void makeBox(SimuData *simu, double *efac, SimuData *&boxed, generateMock_info& args_info)
{
  uint32_t goodParticles = 0;
  double ranges[3][2] = {
    { args_info.rangeX_min_arg, args_info.rangeX_max_arg },
    { args_info.rangeY_min_arg, args_info.rangeY_max_arg },
    { args_info.rangeZ_min_arg, args_info.rangeZ_max_arg }
  };
  double mul[3];
  int *particle_id;

  boxed = new SimuData;
  boxed->Hubble = simu->Hubble;
  boxed->Omega_M = simu->Omega_M;
  boxed->Omega_Lambda = simu->Omega_Lambda;
  boxed->time = simu->time;
  boxed->BoxSize = simu->BoxSize;

  for (uint32_t i = 0; i < simu->NumPart; i++)
    {
      bool acceptance = true;

      for (int j = 0; j < 3; j++)
	acceptance = 
	  acceptance &&
	  (simu->Pos[j][i] > ranges[j][0]) && 
	  (simu->Pos[j][i] < ranges[j][1]);
      
      if (acceptance)
	goodParticles++;
    }
  
  for (int j = 0; j < 3; j++)
    {
      boxed->Pos[j] = new float[goodParticles];
      boxed->Vel[j] = 0;
      mul[j] = 1.0/(ranges[j][1] - ranges[j][0]);
    }

  boxed->NumPart = goodParticles;

  particle_id = new int[goodParticles];
  double *expansion_fac = new double[goodParticles];

  uint32_t k = 0;
  for (uint32_t i = 0; i < simu->NumPart; i++)
    {
      bool acceptance = true;

      for (int j = 0; j < 3; j++)
	acceptance = 
	  acceptance &&
	  (simu->Pos[j][i] > ranges[j][0]) && 
	  (simu->Pos[j][i] < ranges[j][1]);
      
      if (acceptance)
	{
	  for (int j = 0; j < 3; j++)
	    {
	      boxed->Pos[j][k] = (simu->Pos[j][i]-ranges[j][0])*mul[j];
	      assert(boxed->Pos[j][k] > 0);
	      assert(boxed->Pos[j][k] < 1);
	    }
	  particle_id[k] = i;
          expansion_fac[k] = efac[i];
	  k++;
	}
    }

  NcFile f(args_info.outputParameter_arg, NcFile::Replace);

  f.add_att("range_x_min", ranges[0][0]);
  f.add_att("range_x_max", ranges[0][1]);
  f.add_att("range_y_min", ranges[1][0]);
  f.add_att("range_y_max", ranges[1][1]);
  f.add_att("range_z_min", ranges[2][0]);
  f.add_att("range_z_max", ranges[2][1]);

  NcDim *NumPart_dim = f.add_dim("numpart_dim", boxed->NumPart);
  NcVar *v = f.add_var("particle_ids", ncInt, NumPart_dim);
  NcVar *v2 = f.add_var("expansion", ncDouble, NumPart_dim);

  v->put(particle_id, boxed->NumPart);
  v2->put(expansion_fac, boxed->NumPart);

  delete[] particle_id;
  delete[] expansion_fac;
}

int main(int argc, char **argv)
{
  generateMock_info args_info;
  generateMock_conf_params args_params;
  SimuData *simu, *simuOut;
 
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

  if (args_info.ramsesBase_given || args_info.ramsesId_given)
    {
      if (args_info.ramsesBase_given && args_info.ramsesId_given)  
	simu = doLoadRamses(args_info.ramsesBase_arg, 
			    args_info.ramsesId_arg, 
			    args_info.axis_arg, false);
      else
	{
	  cerr << "Both ramsesBase and ramsesId are required to be able to load snapshots" << endl;
	  return 1;
	}

      if (simu == 0)
	{
	  cerr << "Error while loading" << endl;
	  return 1;
	}
    }
  else if (args_info.gadget_given)
    {
      simu = doLoadGadget(args_info.gadget_arg, args_info.axis_arg, false);      
      if (simu == 0)
	{
	  cerr << "Error while loading " << endl;
	  return 1;
	}
    }
  else
    {
      cerr << "Either a ramses snapshot or a gadget snapshot is required." << endl;
      return 1;
    }
  cout << "Hubble = " << simu->Hubble << endl;
  cout << "Boxsize = " << simu->BoxSize << endl;
  cout << "Omega_M = " << simu->Omega_M << endl;
  cout << "Omega_Lambda = " << simu->Omega_Lambda << endl;

  double *expfact;

  if (args_info.cosmo_flag)
    metricTransform(simu, args_info.axis_arg, args_info.preReShift_flag, args_info.peculiarVelocities_flag, expfact);
  else
    { 
      expfact = new double[simu->NumPart];
      for (int j = 0; j < simu->NumPart; j++) expfact[j] = 1.0;
    }

  makeBox(simu, expfact, simuOut, args_info);
  delete simu;

  generateOutput(simuOut, args_info.axis_arg, args_info.output_arg);

  delete simuOut;
  
  return 0;
}
