#include <netcdfcpp.h>
#include <CosmoTool/fortran.hpp>
#include "particleInfo.hpp"

using namespace std;
using namespace CosmoTool;

bool loadParticleInfo(ParticleInfo& info,
		      const std::string& particles, 
		      const std::string& extra_info)
{
  try
    {
      UnformattedRead f(particles);

      int numpart;
      
      f.beginCheckpoint();
      numpart = f.readInt32();
      f.endCheckpoint();
      
      info.particles.resize(numpart);
      
      f.beginCheckpoint();
      for (int i = 0; i < numpart; i++)
	info.particles[i].x = f.readReal32();
      f.endCheckpoint();
      
      f.beginCheckpoint();
      for (int i = 0; i < numpart; i++)
	info.particles[i].y = f.readReal32();
      f.endCheckpoint();
      
      f.beginCheckpoint();
      for (int i = 0; i < numpart; i++)
	info.particles[i].z = f.readReal32();
      f.endCheckpoint();
    }
  catch (const NoSuchFileException& e)
    {
      return false;
    }

  NcFile f_info(extra_info.c_str());
  
  if (!f_info.is_valid())
    return false;

  info.ranges[0][0] = f_info.get_att("range_x_min")->as_double(0);
  info.ranges[0][1] = f_info.get_att("range_x_max")->as_double(0);
  info.ranges[1][0] = f_info.get_att("range_y_min")->as_double(0);
  info.ranges[1][1] = f_info.get_att("range_y_max")->as_double(0);
  info.ranges[2][0] = f_info.get_att("range_z_min")->as_double(0);
  info.ranges[2][1] = f_info.get_att("range_z_max")->as_double(0);

  return true;
}
