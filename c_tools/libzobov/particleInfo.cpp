#include <netcdfcpp.h>
#include <CosmoTool/fortran.hpp>
#include "particleInfo.hpp"

using namespace std;
using namespace CosmoTool;

bool loadParticleInfo(ParticleInfo& info,
		      const std::string& particles, 
		      const std::string& extra_info)
{
      int numpart;
      int isObservation;

  NcFile f_info(extra_info.c_str());
 
  if (!f_info.is_valid())
    return false;

  info.ranges[0][0] = f_info.get_att("range_x_min")->as_double(0);
  info.ranges[0][1] = f_info.get_att("range_x_max")->as_double(0);
  info.ranges[1][0] = f_info.get_att("range_y_min")->as_double(0);
  info.ranges[1][1] = f_info.get_att("range_y_max")->as_double(0);
  info.ranges[2][0] = f_info.get_att("range_z_min")->as_double(0);
  info.ranges[2][1] = f_info.get_att("range_z_max")->as_double(0);
  info.mask_index = f_info.get_att("mask_index")->as_int(0); //PMS
  isObservation = f_info.get_att("is_observation")->as_int(0); //PMS

  for (int i = 0; i < 3; i++)
    info.length[i] = info.ranges[i][1] - info.ranges[i][0];

  try
    {
      UnformattedRead f(particles);

      float mul, offset;
      
      f.beginCheckpoint();
      numpart = f.readInt32();
      f.endCheckpoint();
      
      info.particles.resize(numpart);
      
      offset = info.ranges[0][0];
      mul = info.ranges[0][1] - info.ranges[0][0];
      f.beginCheckpoint();
      for (int i = 0; i < numpart; i++)
	info.particles[i].x = mul*f.readReal32();
      f.endCheckpoint();
      
      offset = info.ranges[1][0];
      mul = info.ranges[1][1] - info.ranges[1][0];
      f.beginCheckpoint();
      for (int i = 0; i < numpart; i++)
	info.particles[i].y = mul*f.readReal32();
      f.endCheckpoint();
      
      offset = info.ranges[2][0];
      mul = info.ranges[2][1] - info.ranges[2][0];
      f.beginCheckpoint();
      for (int i = 0; i < numpart; i++)
	info.particles[i].z = mul*f.readReal32();
      f.endCheckpoint();

      if (!isObservation) {
        for (int i = 0; i < numpart; i++) {
	        info.particles[i].x += info.ranges[0][0];
	        info.particles[i].y += info.ranges[1][0];
	        info.particles[i].z += info.ranges[2][0];
        }
      }
    }
  catch (const NoSuchFileException& e)
    {
      return false;
    }

  return true;
}
