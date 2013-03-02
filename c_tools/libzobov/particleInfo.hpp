/*+
    VIDE -- Void IDEntification pipeline -- ./c_tools/libzobov/particleInfo.hpp
    Copyright (C) 2010-2013 Guilhem Lavaux
    Copyright (C) 2011-2013 Paul M. Sutter

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
+*/

#ifndef _PARTICLE_INFO_HEADER_HPP
#define _PARTICLE_INFO_HEADER_HPP

#include <vector>
#include <string>

struct ParticleData {
  float x, y, z;
};

typedef std::vector<ParticleData> ParticleVector;

struct ParticleInfo
{
  ParticleVector particles;
  float ranges[3][2];
  float length[3];
  int   mask_index; // PMS
};

bool loadParticleInfo(ParticleInfo& info,
		      const std::string& particles, 
		      const std::string& extra_info);

#endif
