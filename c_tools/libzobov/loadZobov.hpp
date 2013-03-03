/*+
    VIDE -- Void IDEntification pipeline -- ./c_tools/libzobov/loadZobov.hpp
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



#ifndef __LOAD_ZOBOV_HPP
#define __LOAD_ZOBOV_HPP

#include <vector>

struct ZobovZone
{
  std::vector<int> pId;
};

struct ZobovVoid
{
  std::vector<int> zId;
  float proba;
  int numParticles, coreParticle;
  float volume;
  float barycenter[3];
  float nearestBoundary;
};

struct ZobovRep
{
  std::vector<ZobovZone> allZones;
  std::vector<ZobovVoid> allVoids;
  std::vector<float> particleVolume;
};

struct ZobovParticle
{
  float x, y, z;
};

bool loadZobov(const char *descName,
	       const char *adjName, const char *voidName,
	       const char *volName, ZobovRep& z);

bool loadZobovParticles(const char *fname, std::vector<ZobovParticle>& particles);

#endif
