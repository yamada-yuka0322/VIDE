#include <cmath>
#include <cassert>
#include <iostream>
#include <fstream>
#include <string>
#include <CosmoTool/loadSimu.hpp>
#include <CosmoTool/loadRamses.hpp>
#include <CosmoTool/loadGadget.hpp>
#include <CosmoTool/loadFlash.hpp>
#include <CosmoTool/interpolate.hpp>
#include <CosmoTool/fortran.hpp>
#include "generateMock_conf.h"
#include "gslIntegrate.hpp"
#include <netcdfcpp.h>

using namespace std;
using namespace CosmoTool;

#define LIGHT_SPEED 299792.458

static double gadgetUnit=1e-3;

SimuData *doLoadRamses(const char *basename, int baseid, int velAxis, bool goRedshift)
{
  SimuData *d, *outd;

  d = loadRamsesSimu(basename, baseid, -1, true, 0);
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
  while (d = loadRamsesSimu(basename, baseid, curCpu, true, NEED_POSITION|NEED_VELOCITY|NEED_GADGET_ID))
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



SimuData *myLoadGadget(const char *fname, int id, int flags)
{
  SimuData *sim = loadGadgetMulti(fname, id, flags);
  sim->BoxSize *= gadgetUnit*1000;
  for (int j = 0; j < 3; j++)
  {
    if (sim->Pos[j] != 0) {
      for (long i = 0; i < sim->NumPart; i++)
        sim->Pos[j][i] *= gadgetUnit*1000;
    }
  }
  return sim;
}

SimuData *doLoadSimulation(const char *gadgetname, int velAxis, bool goRedshift, SimuData *(*loadFunction)(const char *fname, int id, int flags))
{
  SimuData *d, *outd;
  bool singleFile = false;

  try
    {
      d = loadFunction(gadgetname, -1, 0);
      singleFile = true;
    }
  catch (const NoSuchFileException& e)
    {
      try
	{
	  d = loadFunction(gadgetname, 0, 0);
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

  int curCpu = singleFile ? -1 : 0;
  cout << "loading file 0 " << endl;
  try
    {
      while (1)
	{
	  d = loadFunction(gadgetname, curCpu, NEED_POSITION|NEED_VELOCITY|NEED_GADGET_ID);
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
          if (singleFile)
            break;
	  curCpu++;
	  cout << "loading file " << curCpu  << endl;
	}
    }
  catch (const NoSuchFileException& e)
    {
    }

  return outd;
}


SimuData *doLoadMultidark(const char *multidarkname)
{
  SimuData *outd;
  FILE *fp;
  int actualNumPart;

  outd = new SimuData;
  cout << "opening multidark file " << multidarkname << endl;
  fp = fopen(multidarkname, "r");
  if (fp == NULL) {
    cout << "could not open file!" << endl;
    return 0;
  }
  fscanf(fp, "%f\n", &outd->BoxSize);
  fscanf(fp, "%f\n", &outd->Omega_M);
  fscanf(fp, "%f\n", &outd->Hubble);
  fscanf(fp, "%f\n", &outd->time);
  fscanf(fp, "%ld\n", &outd->NumPart);

  outd->time = 1./(1.+outd->time); // convert to scale factor
  outd->TotalNumPart = outd->NumPart;
  outd->Omega_Lambda = 1.0 - outd->Omega_M;

  for (int k = 0; k < 3; k++)
    outd->Pos[k] = new float[outd->NumPart];
  outd->Vel[2] = new float[outd->NumPart];
  outd->Id = new int[outd->NumPart];

  cout << "loading multidark particles" << endl;
  actualNumPart = 0;
	for (int i = 0; i < outd->NumPart; i++) {
    fscanf(fp, "%d %d %f %f %f\n", &outd->Id[i], 
                                &outd->Pos[0][i], &outd->Pos[1][i], 
                                &outd->Pos[2][i], &outd->Vel[2][i]);

    if (outd->Id[i] == -99 && 
        outd->Pos[0][i] == -99 && outd->Pos[1][i] == -99 && 
        outd->Pos[2][i] == -99 && outd->Vel[2][i] == -99) {
      break;
    } else {
      actualNumPart++;
    }
  }
  fclose(fp);

  outd->NumPart = actualNumPart;
  outd->TotalNumPart = actualNumPart;
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

Interpolate make_cosmological_redshift(double OM, double OL, double z0, double z1, int N = 5000)
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

void metricTransform(SimuData *data, int axis, bool reshift, bool pecvel, double*& expfact, bool cosmo_flag)
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

  double z0 = 1/data->time - 1;
  Interpolate z_vs_D = 
    make_cosmological_redshift(data->Omega_M, data->Omega_Lambda, 
			       0., z0+8*data->BoxSize*100/LIGHT_SPEED);
  // Redshift 2*z0 should be sufficient ? This is fragile. 
  //A proper solver is needed here.
  double z_base = reshift ? z0 : 0;
  TotalExpansion e_computer;
  double baseComovingDistance;

  expfact = new double[data->NumPart];

  cout << "Using base redshift z=" << z0 << " " << z0+8*data->BoxSize*100/LIGHT_SPEED << endl;

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
          else if (cosmo_flag)
	    z = (z_vs_D.compute(reduced_red)-z_base)*LIGHT_SPEED/100.;
          else
            z = reduced_red*LIGHT_SPEED/100.0;

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

  cout << "Writing RA..." << endl;
  f.beginCheckpoint();
  for (uint32_t i = 0; i < data->NumPart; i++)
    {
      f.writeReal32(data->Id[i]);
    }
  f.endCheckpoint();

  cout << "Writing Dec..." << endl;
  f.beginCheckpoint();
  for (uint32_t i = 0; i < data->NumPart; i++)
    {
      f.writeReal32(data->Id[i]);
    }
  f.endCheckpoint();

  cout << "Writing redshift..." << endl;
  f.beginCheckpoint();
  for (uint32_t i = 0; i < data->NumPart; i++)
    {
      f.writeReal32(data->Id[i]);
    }
  f.endCheckpoint();

  cout << "Writing unique ID..." << endl;
  f.beginCheckpoint();
  for (uint32_t i = 0; i < data->NumPart; i++)
    {
      f.writeReal32(data->Id[i]);
    }
  f.endCheckpoint();


}

void makeBox(SimuData *simu, double *efac, SimuData *&boxed, generateMock_info& args_info)
{
  float subsample = args_info.subsample_given ? args_info.subsample_arg : 1.0;
  uint32_t goodParticles = 0;
  double ranges[3][2] = {
    { args_info.rangeX_min_arg, args_info.rangeX_max_arg },
    { args_info.rangeY_min_arg, args_info.rangeY_max_arg },
    { args_info.rangeZ_min_arg, args_info.rangeZ_max_arg }
  };
  double mul[3];
  float minmax[2][3];
  int *particle_id;
  bool *random_acceptance = 0;

  boxed = new SimuData;
  boxed->Hubble = simu->Hubble;
  boxed->Omega_M = simu->Omega_M;
  boxed->Omega_Lambda = simu->Omega_Lambda;
  boxed->time = simu->time;
  boxed->BoxSize = simu->BoxSize;

  random_acceptance = new bool[simu->NumPart];

  for (int j = 0; j < 3; j++) minmax[1][j] = minmax[0][j] = simu->Pos[j][0];

  for (uint32_t i = 0; i < simu->NumPart; i++)
    {
      bool acceptance = true;

      for (int j = 0; j < 3; j++) {
	acceptance = 
	  acceptance &&
	  (simu->Pos[j][i] > ranges[j][0]) && 
	  (simu->Pos[j][i] < ranges[j][1]);
        minmax[0][j] = min(simu->Pos[j][i], minmax[0][j]);
        minmax[1][j] = max(simu->Pos[j][i], minmax[1][j]);
      }
      random_acceptance[i] = acceptance && (drand48() <= subsample);
      if (random_acceptance[i])
	goodParticles++;
    }

  cout << "Subsample fraction: " << subsample << endl;
  cout << "Min range = " << ranges[0][0] << " " << ranges[1][0] << " " << ranges[2][0] << endl;
  cout << "Max range = " << ranges[0][1] << " " << ranges[1][1] << " " << ranges[2][1] << endl;

  cout << "Min position = " << minmax[0][0] << " " << minmax[0][1] << " " << minmax[0][2] << endl;
  cout << "Max position = " << minmax[1][0] << " " << minmax[1][1] << " " << minmax[1][2] << endl;
 
  cout << "Number of accepted particles: " << goodParticles << endl; 

  for (int j = 0; j < 3; j++)
    {
      boxed->Pos[j] = new float[goodParticles];
      boxed->Vel[j] = 0;
      mul[j] = 1.0/(ranges[j][1] - ranges[j][0]);
    }

  cout << "Rescaling factors = " << mul[0] << " " << mul[1] << " " << mul[2] << endl;
  boxed->NumPart = goodParticles;

  particle_id = new int[goodParticles];
  double *expansion_fac = new double[goodParticles];

  uint32_t k = 0;
  for (uint32_t i = 0; i < simu->NumPart; i++)
    {
      bool acceptance = random_acceptance[i]; 
      
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

  delete[] random_acceptance;

  NcFile f(args_info.outputParameter_arg, NcFile::Replace);

  f.add_att("range_x_min", ranges[0][0]);
  f.add_att("range_x_max", ranges[0][1]);
  f.add_att("range_y_min", ranges[1][0]);
  f.add_att("range_y_max", ranges[1][1]);
  f.add_att("range_z_min", ranges[2][0]);
  f.add_att("range_z_max", ranges[2][1]);
  f.add_att("mask_index", -1);

  NcDim *NumPart_dim = f.add_dim("numpart_dim", boxed->NumPart);
  NcVar *v = f.add_var("particle_ids", ncInt, NumPart_dim);
  NcVar *v2 = f.add_var("expansion", ncDouble, NumPart_dim);

  v->put(particle_id, boxed->NumPart);
  v2->put(expansion_fac, boxed->NumPart);

  delete[] particle_id;
  delete[] expansion_fac;


  FILE *fp = fopen("sample_info.txt", "w");
  fprintf(fp, "x_min = %f\n", ranges[0][0]);
  fprintf(fp, "x_max = %f\n", ranges[0][1]);
  fprintf(fp, "y_min = %f\n", ranges[1][0]);
  fprintf(fp, "y_max = %f\n", ranges[1][1]);
  fprintf(fp, "z_min = %f\n", ranges[2][0]);
  fprintf(fp, "z_max = %f\n", ranges[2][1]);
  fprintf(fp, "mask_index = -1\n");
  fprintf(fp, "total_particles = %d\n", boxed->NumPart);
  fclose(fp);
}

void makeBoxFromParameter(SimuData *simu, double *efac, SimuData* &boxed, generateMock_info& args_info)
{
  NcFile f(args_info.inputParameter_arg);
  NcVar *v;
  int *particle_id;
  double *expansion_fac;

  boxed = new SimuData;
  boxed->Hubble = simu->Hubble;
  boxed->Omega_M = simu->Omega_M;
  boxed->Omega_Lambda = simu->Omega_Lambda;
  boxed->time = simu->time;
  boxed->BoxSize = simu->BoxSize;

  NcVar *v_id = f.get_var("particle_ids");
  long *edges1;
  double ranges[3][2];
  double mul[3];

  edges1 = v_id->edges();
  assert(v_id->num_dims()==1);

  boxed->NumPart = edges1[0];
  delete[] edges1;

  particle_id = new int[boxed->NumPart];

  v_id->get(particle_id, boxed->NumPart);

  ranges[0][0] = f.get_att("range_x_min")->as_double(0);
  ranges[0][1] = f.get_att("range_x_max")->as_double(0);
  ranges[1][0] = f.get_att("range_y_min")->as_double(0);
  ranges[1][1] = f.get_att("range_y_max")->as_double(0);
  ranges[2][0] = f.get_att("range_z_min")->as_double(0);
  ranges[2][1] = f.get_att("range_z_max")->as_double(0);

  for (int j = 0; j < 3; j++)
    {
      boxed->Pos[j] = new float[boxed->NumPart];
      boxed->Vel[j] = 0;
      mul[j] = 1.0/(ranges[j][1] - ranges[j][0]);
    }
  
  uint32_t k = 0;
  for (uint32_t i = 0; i < boxed->NumPart; i++)
    {
      int id = particle_id[i];

      for (int j = 0; j < 3; j++)
	{
	  boxed->Pos[j][i] = (simu->Pos[j][id]-ranges[j][0])*mul[j];
	}
    } 

  delete[] particle_id;
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

  gadgetUnit=args_info.gadgetUnit_arg;

  if (args_info.ramsesBase_given || args_info.ramsesId_given)
    {
      if (args_info.ramsesBase_given && args_info.ramsesId_given)   {
	simu = doLoadRamses(args_info.ramsesBase_arg, 
			    args_info.ramsesId_arg, 
			    args_info.axis_arg, false);
      }
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
  else if (args_info.gadget_given || args_info.flash_given || args_info.multidark_given)
    {
      if (args_info.gadget_given && args_info.flash_given)
	{
	  cerr << "Do not know which file to use: Gadget or Flash ?" << endl;
	  return 1;
	}

      if (args_info.multidark_given) {
        simu = doLoadMultidark(args_info.multidark_arg);      
      }

      if (args_info.gadget_given) {
	simu = doLoadSimulation(args_info.gadget_arg, args_info.axis_arg, false, myLoadGadget);      
      }
      if (args_info.flash_given) {
	simu = doLoadSimulation(args_info.flash_arg, args_info.axis_arg, false, loadFlashMulti); 
      }

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

  metricTransform(simu, args_info.axis_arg, args_info.preReShift_flag,
		  args_info.peculiarVelocities_flag, expfact, 
		  args_info.cosmo_flag);

  if (args_info.inputParameter_given)
    makeBoxFromParameter(simu, expfact, simuOut, args_info);
  else
    makeBox(simu, expfact, simuOut, args_info);

  delete simu;

  generateOutput(simuOut, args_info.axis_arg, args_info.output_arg);

  delete simuOut;
 
  printf("Done!\n"); 
  return 0;
}
