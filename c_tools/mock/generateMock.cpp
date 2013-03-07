/*+
    VIDE -- Void IDEntification pipeline -- ./c_tools/mock/generateMock.cpp
    Copyright (C) 2010-2013 Guilhem Lavaux
    Copyright (C) 2011-2013 P. M. Sutter

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



#include <gsl/gsl_rng.h>
#include <boost/function.hpp>
#include <boost/bind.hpp>
#include <boost/format.hpp>
#include <cmath>
#include <cassert>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>
#include <CosmoTool/loadSimu.hpp>
#include <CosmoTool/interpolate.hpp>
#include <CosmoTool/fortran.hpp>
#include <CosmoTool/algo.hpp>
#include "generateMock_conf.h"
#include "gslIntegrate.hpp"
#include <netcdfcpp.h>
#include "simulation_loader.hpp"

using namespace std;
using namespace CosmoTool;
using boost::format;

#define LIGHT_SPEED 299792.458

typedef boost::function2<void, SimuData*, double*> MetricFunctor;

struct TotalExpansion
{
  double Omega_M, Omega_L;

  double operator()(double z)
  {
    return 1/sqrt(Omega_M*cube(1+z) + Omega_L);
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

void metricTransform(SimuData *data, int axis, bool reshift, bool pecvel, double* expfact, bool cosmo_flag)
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

          if (expfact)
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
      f.writeReal32(data->Pos[x0][i]);
    }
  f.endCheckpoint();

  cout << "Writing Dec..." << endl;
  f.beginCheckpoint();
  for (uint32_t i = 0; i < data->NumPart; i++)
    {
      f.writeReal32(data->Pos[x1][i]);
    }
  f.endCheckpoint();

  cout << "Writing redshift..." << endl;
  f.beginCheckpoint();
  for (uint32_t i = 0; i < data->NumPart; i++)
    {
      f.writeReal32(data->Pos[x2][i]*LIGHT_SPEED);
    }
  f.endCheckpoint();

  long *uniqueID = data->as<long>("uniqueID");
  if (uniqueID != 0)
    {
      cout << "Writing unique ID..." << endl;
      f.beginCheckpoint();
      for (uint32_t i = 0; i < data->NumPart; i++)
	{
	  f.writeInt64(uniqueID[i]);
	}
      f.endCheckpoint();
    }
}

// This function prepares the list of targets for the specified snapshot. The target list is not
// cleared. That way new particles can be appended if this is a multi-file snapshot.
void selectBox(SimuData *simu, std::vector<long>& targets, generateMock_info& args_info, SimulationPreprocessor *preselect)
{
  double ranges[3][2] = {
    { args_info.rangeX_min_arg, args_info.rangeX_max_arg },
    { args_info.rangeY_min_arg, args_info.rangeY_max_arg },
    { args_info.rangeZ_min_arg, args_info.rangeZ_max_arg }
  };
  long numAccepted =0;

  for (uint32_t i = 0; i < simu->NumPart; i++)
    {
      bool acceptance = true;
      SingleParticle p;

      for (int j = 0; j < 3; j++) {
	acceptance = 
	  acceptance &&
	  (simu->Pos[j][i] > ranges[j][0]) && 
	  (simu->Pos[j][i] <= ranges[j][1]);	
        p.Pos[j] = simu->Pos[j][i];
        p.Vel[j] = (simu->Vel[j] != 0) ? simu->Vel[j][i] : 0;
      }
      p.ID = (simu->Id != 0) ? simu->Id[i] : -1;

      if (preselect != 0)
        acceptance = acceptance && preselect->accept(p);
      
      if (acceptance) {
	targets.push_back(i);
        numAccepted++; 
      }
    }
  cout << "Num accepted here = " << numAccepted <<  " / input = " << simu->NumPart << endl; 
}

class PreselectParticles: public SimulationPreprocessor
{
private:
  gsl_rng *rng;
  double subsample;
  int seed;
public:
  PreselectParticles(double s, int seed_value)
    : subsample(s), seed(seed_value), rng(gsl_rng_alloc(gsl_rng_default))
  {
    gsl_rng_set(rng, seed);
  }

  virtual ~PreselectParticles()
  {
    gsl_rng_free(rng);
  } 

  bool accept(const SingleParticle& p)
  {
    return gsl_rng_uniform(rng) < subsample;
  }

  void reset()
  {
    gsl_rng_set(rng, seed);
  }
};

void createBox(SimuData *simu, vector<long>& targets, vector<long>& snapshot_split, SimuData *& boxed, generateMock_info& args_info)
{
  double *ranges = new double[6];
  double *mul = new double[3];
  long *simu_uniqueID = simu->as<long>("uniqueID");

  ranges[0] = args_info.rangeX_min_arg;
  ranges[1] = args_info.rangeX_max_arg;
  ranges[2] = args_info.rangeY_min_arg;
  ranges[3] = args_info.rangeY_max_arg;
  ranges[4] = args_info.rangeZ_min_arg;
  ranges[5] = args_info.rangeZ_max_arg;

  boxed = new SimuData;
  boxed->Hubble = simu->Hubble;
  boxed->Omega_M = simu->Omega_M;
  boxed->Omega_Lambda = simu->Omega_Lambda;
  boxed->time = simu->time;
  boxed->BoxSize = simu->BoxSize;
  boxed->NumPart = targets.size();

  for (int j = 0; j < 3; j++)
    {
      boxed->Pos[j] = new float[boxed->NumPart];
      boxed->Vel[j] = 0;
      mul[j] = 1.0/(ranges[2*j+1] - ranges[2*j+0]);
    }
  cout << "Min range = " << ranges[0] << " " << ranges[2] << " " << ranges[4] << endl;
  cout << "Max range = " << ranges[1] << " " << ranges[3] << " " << ranges[5] << endl;
  cout << "Number of accepted particles: " << boxed->NumPart << endl; 
  cout << "Rescaling factors = " << mul[0] << " " << mul[1] << " " << mul[2] << endl;

  long *uniqueID = new long[boxed->NumPart];
  long *particle_id = new long[boxed->NumPart];
  double *expansion_fac = new double[boxed->NumPart];
  long *snap_split = new long[snapshot_split.size()];
  int *numsnap_info = new int[1];

  copy(targets.begin(), targets.end(), particle_id);
  copy(snapshot_split.begin(), snapshot_split.end(), snap_split);
  *numsnap_info = snapshot_split.size();

  boxed->new_attribute("particle_id", particle_id, delete_adaptor<long>);
  boxed->new_attribute("expansion_fac", expansion_fac, delete_adaptor<double>);
  boxed->new_attribute("uniqueID", uniqueID, delete_adaptor<long>); 
  boxed->new_attribute("mul", mul, delete_adaptor<double>);
  boxed->new_attribute("ranges", ranges, delete_adaptor<double>);
  boxed->new_attribute("snapshot_split", snap_split, delete_adaptor<long>);
  boxed->new_attribute("num_snapshots", numsnap_info, delete_adaptor<int>);
}

void buildBox(SimuData *simu, long num_targets, long loaded, 
	      SimuData *boxed, double *efac)
{
  uint32_t k = 0;
  long *uniqueID = boxed->as<long>("uniqueID");
  long *simu_uniqueID = simu->as<long>("uniqueID");
  double *expansion_fac = boxed->as<double>("expansion_fac");
  double *mul = boxed->as<double>("mul");
  double *ranges = boxed->as<double>("ranges");
  long *particle_id = boxed->as<long>("particle_id");

  for (uint32_t i = 0; i < num_targets; i++, loaded++)
    {
      long pid = particle_id[loaded];
      //assert(pid < simu->NumPart);
      assert(loaded < boxed->NumPart);
      
      for (int j = 0; j < 3; j++)
	{
	  boxed->Pos[j][loaded] = max(min((simu->Pos[j][pid]-ranges[j*2])*mul[j], double(1)), double(0));
	  assert(boxed->Pos[j][loaded] >= 0);
	  assert(boxed->Pos[j][loaded] <= 1);
	}
      uniqueID[loaded] = (simu_uniqueID != 0) ? simu_uniqueID[pid] : 0;
      expansion_fac[loaded] = efac[pid];
    }
}

void saveBox(SimuData *&boxed, const std::string& outbox)
{
  double *ranges = boxed->as<double>("ranges");
  NcFile f(outbox.c_str(), NcFile::Replace, 0, 0, NcFile::Netcdf4);
  long *particle_id = boxed->as<long>("particle_id");
  double *expansion_fac = boxed->as<double>("expansion_fac");
  long *snapshot_split = boxed->as<long>("snapshot_split");
  int num_snapshots = *boxed->as<int>("num_snapshots");
  long *uniqueID = boxed->as<long>("uniqueID");

  f.add_att("range_x_min", ranges[0]);
  f.add_att("range_x_max", ranges[1]);
  f.add_att("range_y_min", ranges[2]);
  f.add_att("range_y_max", ranges[3]);
  f.add_att("range_z_min", ranges[4]);
  f.add_att("range_z_max", ranges[5]);
  f.add_att("mask_index", -1);
  f.add_att("is_observation", 0);

  NcDim *NumPart_dim = f.add_dim("numpart_dim", boxed->NumPart);
  NcDim *NumSnap_dim = f.add_dim("numsnap_dim", num_snapshots);
  NcVar *v = f.add_var("particle_ids", ncLong, NumPart_dim);
  NcVar *v2 = f.add_var("expansion", ncDouble, NumPart_dim);
  NcVar *v3 = f.add_var("snapshot_split", ncLong, NumSnap_dim);

  v->put(particle_id, boxed->NumPart);
  v2->put(expansion_fac, boxed->NumPart);  
  v3->put(snapshot_split, num_snapshots);
  if (uniqueID != 0)
    {
      NcVar *v4 = f.add_var("unique_ids_lsb", ncLong, NumPart_dim);
      NcVar *v5 = f.add_var("unique_ids_msb", ncLong, NumPart_dim);

      nclong *tmp_int = new nclong[boxed->NumPart];
      for (long i = 0; i < boxed->NumPart; i++)
         tmp_int[i] = (nclong)(((unsigned long)uniqueID[i]) & 0xffffffff);
      v4->put(tmp_int, boxed->NumPart);
      for (long i = 0; i < boxed->NumPart; i++)
         tmp_int[i] = (nclong)((((unsigned long)uniqueID[i]) & 0xffffffff) >> 32);
      v5->put(tmp_int, boxed->NumPart);
      delete[] tmp_int;
    }
}

void makeBoxFromParameter(SimuData *simu, SimuData* &boxed, generateMock_info& args_info)
{
  NcFile f(args_info.inputParameter_arg);
  NcVar *v;
  long *particle_id;
  double *expansion_fac;
  long *uniqueID;
  long *snapshot_split;
  int *num_snapshots = new int[1];

  boxed = new SimuData;
  boxed->Hubble = simu->Hubble;
  boxed->Omega_M = simu->Omega_M;
  boxed->Omega_Lambda = simu->Omega_Lambda;
  boxed->time = simu->time;
  boxed->BoxSize = simu->BoxSize;

  NcVar *v_id = f.get_var("particle_ids");
  NcVar *v_snap = f.get_var("snapshot_split");
  long *edges1;
  long *dim_snap;
  double *ranges;
  double *mul;

  edges1 = v_id->edges();
  dim_snap = v_snap->edges();
  assert(v_id->num_dims()==1);
  assert(v_snap->num_dims()==1);

  boxed->NumPart = edges1[0];
  *num_snapshots = dim_snap[0];
  delete[] dim_snap;
  delete[] edges1;

  particle_id = new long[boxed->NumPart];
  uniqueID = new long[boxed->NumPart];
  mul = new double[3];
  ranges = new double[6];
  snapshot_split = new long[*num_snapshots];
  expansion_fac = new double[boxed->NumPart];


  boxed->new_attribute("uniqueID", uniqueID, delete_adaptor<long>);
  boxed->new_attribute("mul", mul, delete_adaptor<double>);
  boxed->new_attribute("ranges", ranges, delete_adaptor<double>);
  boxed->new_attribute("particle_id", particle_id, delete_adaptor<long>);
  boxed->new_attribute("num_snapshots", num_snapshots, delete_adaptor<int>);
  boxed->new_attribute("snapshot_split", snapshot_split, delete_adaptor<long>);
  boxed->new_attribute("expansion_fac", expansion_fac, delete_adaptor<double>);

  v_id->get(particle_id, boxed->NumPart);
  v_snap->get(snapshot_split, *num_snapshots);

  ranges[0] = f.get_att("range_x_min")->as_double(0);
  ranges[1] = f.get_att("range_x_max")->as_double(0);
  ranges[2] = f.get_att("range_y_min")->as_double(0);
  ranges[3] = f.get_att("range_y_max")->as_double(0);
  ranges[4] = f.get_att("range_z_min")->as_double(0);
  ranges[5] = f.get_att("range_z_max")->as_double(0);

  for (int j = 0; j < 3; j++)
    {
      boxed->Pos[j] = new float[boxed->NumPart];
      boxed->Vel[j] = 0;
      mul[j] = 1.0/(ranges[2*j+1] - ranges[2*j+0]);
    }
  
  uint32_t k = 0;
  NcVar *v_uniq_lsb = f.get_var("unique_ids_lsb");
  NcVar *v_uniq_msb = f.get_var("unique_ids_lsb");
  nclong *tmp_int;

  tmp_int = new nclong[boxed->NumPart];

  v_uniq_lsb->get(tmp_int, boxed->NumPart);
  for (long i = 0; i < boxed->NumPart; i++)
    uniqueID[i] = tmp_int[i];

  v_uniq_msb->get(tmp_int, boxed->NumPart);
  for (long i = 0; i < boxed->NumPart; i++)
    uniqueID[i] |= (unsigned long)(tmp_int[i]) << 32;

  delete[] tmp_int;

  PreselectParticles *preselect = 0;

  if (args_info.resubsample_given)
    {
      preselect = new PreselectParticles(args_info.resubsample_arg, args_info.resubsample_seed_arg);
      preselect->reset();
    }

  if (preselect == 0)
    return;

  long pid_read = 0, pid_write = 0;
  for (int s_id = 0; s_id < *num_snapshots; s_id++)
    {
      long previous_write = pid_write;
      for (long q = 0; q < snapshot_split[s_id]; q++)
        {
          SingleParticle p;
          p.ID = -1;
          assert(pid_read < boxed->NumPart);
          if (preselect->accept(p))
           {
             particle_id[pid_write] = particle_id[pid_read];
             uniqueID[pid_write] = uniqueID[pid_read];
             expansion_fac[pid_write] = expansion_fac[pid_read];
             pid_write++; 
           }
          pid_read++;
        }
      snapshot_split[s_id] = pid_write - previous_write; 
    }  
  boxed->NumPart = pid_write; 
  delete preselect;
  cout << "Num accepted here = " << pid_write <<  " / input = " << pid_read << endl; 
}

void makeBoxFromSimulation(SimulationLoader *loader, SimuData* &boxed, MetricFunctor metric, generateMock_info& args_info)
{
  vector<long> targets, split;
  long previous_target_num = 0;
  PreselectParticles *preselect = 0;

  if (args_info.resubsample_given)
    {
      preselect = new PreselectParticles(args_info.resubsample_arg, args_info.resubsample_seed_arg);
      preselect->reset();
    }

  for (int nf = 0; nf < loader->num_files(); nf++)
    {
      SimuData *simu;
      double *expfact;

      cout << format("Analyzing and selecting targets in file number %d / %d") % (nf+1) % loader->num_files() << endl;
      simu = loader->loadFile(nf);      

      metric(simu, 0);

      selectBox(simu, targets, args_info, preselect);
      split.push_back(targets.size() - previous_target_num);
      previous_target_num = targets.size();
      
      delete simu;
    }

  createBox(loader->getHeader(), targets, split, boxed, args_info);

  if (preselect)
    delete preselect;
}

int main(int argc, char **argv)
{
  generateMock_info args_info;
  generateMock_conf_params args_params;
  SimuData *simu, *simuOut;
  SimulationLoader *loader;
 
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

  SimulationPreprocessor *preselector = new PreselectParticles(args_info.subsample_arg, args_info.subsample_seed_arg);

  if (args_info.ramsesBase_given || args_info.ramsesId_given)
    {
      if (args_info.ramsesBase_given && args_info.ramsesId_given)   {
	loader = ramsesLoader(args_info.ramsesBase_arg, 
			      args_info.ramsesId_arg, 
			      false,
			      NEED_POSITION|NEED_VELOCITY|NEED_GADGET_ID, preselector);
      }
      else
	{
	  cerr << "Both ramsesBase and ramsesId are required to be able to load snapshots" << endl;
	  return 1;
	}
    }
  else if (args_info.gadget_given)
    {
      loader = gadgetLoader(args_info.gadget_arg, 1/args_info.gadgetUnit_arg, NEED_POSITION|NEED_VELOCITY|NEED_GADGET_ID, preselector);
    }
  else if (args_info.flash_given)
    {
      loader = flashLoader(args_info.flash_arg, NEED_POSITION|NEED_VELOCITY|NEED_GADGET_ID, preselector); 
  }
  else if (args_info.multidark_given)
    {
      loader = multidarkLoader(args_info.multidark_arg, preselector);      
    }
#ifdef SDF_SUPPORT
  else if (args_info.sdf_given)
    {
      loader = sdfLoader(args_info.sdf_arg, NEED_POSITION|NEED_VELOCITY|NEED_GADGET_ID, args_info.sdf_splitting_arg, preselector);
    }
#endif
  else
    {
      cerr << "A simulation snapshot is required to generate a mock catalog." << endl;
      return 1;
    }
     
  if (loader == 0)
    {
      cerr << "Error while loading " << endl;
      return 1;
    }
  simu = loader->getHeader();

  {
    SimuData *header = loader->getHeader();
    cout << "Hubble = " << header->Hubble << endl;
    cout << "Boxsize = " << header->BoxSize << endl;
    cout << "Omega_M = " << header->Omega_M << endl;
    cout << "Omega_Lambda = " << header->Omega_Lambda << endl;
    cout << "Subsample fraction: " << (args_info.subsample_given ? args_info.subsample_arg : 1.0) << endl;
  }
  double *expfact;

  boost::function2<void, SimuData*, double*> metricOperation=
    boost::bind(metricTransform, _1, args_info.axis_arg, args_info.preReShift_flag,
                      args_info.peculiarVelocities_flag, _2, 
                      args_info.cosmo_flag);

  if (args_info.inputParameter_given)
    makeBoxFromParameter(loader->getHeader(), simuOut, args_info);
  else
    makeBoxFromSimulation(loader, simuOut, metricOperation, args_info);

  // Reset the random number generator
  preselector->reset();

  long loaded = 0;
  for (int nf = 0; nf < loader->num_files(); nf++)
    {
      long num_targets = simuOut->as<long>("snapshot_split")[nf];
      cout << format("Building box from particles in %d / %d") % (nf+1) % loader->num_files() << endl;
     
      if (num_targets == 0)
        {
          cout << "No particles selected there. Skipping." << endl;
          continue;
        }
      
      SimuData *simu = loader->loadFile(nf);
      double *efac = new double[simu->NumPart];
      metricOperation(simu, efac);
      buildBox(simu, num_targets, loaded, simuOut, efac);
      
      loaded += num_targets;
      assert(loaded <= simuOut->NumPart);

      delete simu;

      delete[] efac;
    }

  saveBox(simuOut, args_info.outputParameter_arg);
  generateOutput(simuOut, args_info.axis_arg, 
                 args_info.output_arg);
  delete preselector;
 
  double subsample = 1.0;
  if (args_info.subsample_given) subsample = args_info.subsample_arg;
  if (args_info.resubsample_given) subsample = args_info.resubsample_arg;
  printf("Done! %5.2e\n", subsample); 
  return 0;
}
