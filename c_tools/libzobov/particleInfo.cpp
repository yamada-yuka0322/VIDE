/*+
    VIDE -- Void IDEntification pipeline -- ./c_tools/libzobov/particleInfo.cpp
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



#include <cstdlib>
#include <netcdfcpp.h>
#include <CosmoTool/fortran.hpp>
#include "particleInfo.hpp"

using namespace std;
using namespace CosmoTool;

template<bool failure>
double retrieve_attr_safe_double(NcFile& f, const char *name, double defval)
{
  NcAtt *a = f.get_att(name);
  if (a == 0)
    {
      if (failure)
        abort();
      return defval;
    }
  return a->as_double(0);
}

template<bool failure>
int retrieve_attr_safe_int(NcFile& f, const char *name, int defval)
{
  NcAtt *a = f.get_att(name);
  if (a == 0)
    {
      if (failure)
        abort();
      return defval;
    }
  return a->as_int(0);
}

  

bool loadParticleInfo(ParticleInfo& info,
		      const std::string& particles, 
		      const std::string& extra_info)
{
      int numpart;
      int isObservation;

  NcFile f_info(extra_info.c_str());
  NcError nerr(NcError::verbose_nonfatal);
 
  if (!f_info.is_valid())
    return false;

  info.ranges[0][0] = retrieve_attr_safe_double<true>(f_info, "range_x_min", 0);
  info.ranges[0][1] = retrieve_attr_safe_double<true>(f_info, "range_x_max", 0);
  info.ranges[1][0] = retrieve_attr_safe_double<true>(f_info, "range_y_min", 0);
  info.ranges[1][1] = retrieve_attr_safe_double<true>(f_info, "range_y_max", 0);
  info.ranges[2][0] = retrieve_attr_safe_double<true>(f_info, "range_z_min", 0);
  info.ranges[2][1] = retrieve_attr_safe_double<true>(f_info, "range_z_max", 0);
  info.mask_index = retrieve_attr_safe_int<true>(f_info, "mask_index", 0);
  isObservation = retrieve_attr_safe_int<false>(f_info, "is_observation", 0);

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
