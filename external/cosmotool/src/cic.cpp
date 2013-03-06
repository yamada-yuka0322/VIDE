/*+
This is CosmoTool (./src/cic.cpp) -- Copyright (C) Guilhem Lavaux (2007-2013)

guilhem.lavaux@gmail.com

This software is a computer program whose purpose is to provide a toolbox for cosmological
data analysis (e.g. filters, generalized Fourier transforms, power spectra, ...)

This software is governed by the CeCILL license under French law and
abiding by the rules of distribution of free software.  You can  use, 
modify and/ or redistribute the software under the terms of the CeCILL
license as circulated by CEA, CNRS and INRIA at the following URL
"http://www.cecill.info". 

As a counterpart to the access to the source code and  rights to copy,
modify and redistribute granted by the license, users are provided only
with a limited warranty  and the software's author,  the holder of the
economic rights,  and the successive licensors  have only  limited
liability. 

In this respect, the user's attention is drawn to the risks associated
with loading,  using,  modifying and/or developing or reproducing the
software by the user in light of its specific status of free software,
that may mean  that it is complicated to manipulate,  and  that  also
therefore means  that it is reserved for developers  and  experienced
professionals having in-depth computer knowledge. Users are therefore
encouraged to load and test the software's suitability as regards their
requirements in conditions enabling the security of their systems and/or 
data to be ensured and,  more generally, to use and operate it in the 
same conditions as regards security. 

The fact that you are presently reading this means that you have had
knowledge of the CeCILL license and that you accept its terms.
+*/

#include <assert.h>
#include <math.h>
#include <inttypes.h>
#include "cic.hpp"

CICFilter::CICFilter(uint32_t N, double len)
{
  spatialLen = len;
  szGrid = N;
  totalSize = N*N*N;
  densityGrid = new CICType[totalSize];
  resetMesh();
}

CICFilter::~CICFilter()
{
  delete[] densityGrid;
}
  
void CICFilter::resetMesh()
{
  for (uint32_t i = 0; i < totalSize; i++)
    densityGrid[i] = 0;
}

void CICFilter::putParticles(CICParticles *particles, uint32_t N)
{
#if 0
  uint32_t numCorners = 1 << NUMDIMS;

  for (uint32_t i = 0; i < N; i++)
    {
      Coordinates xyz;
      int32_t ixyz[NUMDIMS];
      int32_t rxyz[NUMDIMS];
      CICType alpha[NUMDIMS];
      CICType beta[NUMDIMS];
      for (int j = 0; j < NUMDIMS; j++)
	{
	  xyz[j] = (particles[i].coords[j] / spatialLen * szGrid);
	  ixyz[j] = (int32_t)floor(xyz[j] - 0.5);
	  beta[j] = xyz[j] - ixyz[j] - 0.5;
	  alpha[j] = 1 - beta[j];
	  if (ixyz[j] < 0)
	    ixyz[j] = szGrid-1;
	}

      CICType tot_mass = 0;
      for (int j = 0; j < numCorners; j++)
	{
	  CICType rel_mass = 1;
	  uint32_t idx = 0;
	  uint32_t mul = 1;
	  uint32_t mul2 = 1;

	  for (int k = 0; k < NUMDIMS; k++)
	    {
	      uint32_t ipos = ((j & mul2) != 0);

	      if (ipos == 1)
		{
		  rel_mass *= beta[k];
		}
	      else
		{
		  rel_mass *= alpha[k];
		}

	      rxyz[k] = ixyz[k] + ipos;
	      
	      if (rxyz[k] >= szGrid)
		idx += (rxyz[k] - szGrid) * mul;
	      else
		idx += rxyz[k] * mul;

	      mul2 *= 2;
	      mul *= szGrid;
	    }

	  assert(rel_mass > 0);
	  assert(rel_mass < 1);
	  assert(idx < totalSize);
	  densityGrid[idx] += rel_mass * particles[i].mass;
	  tot_mass += rel_mass;
	}
      assert(tot_mass < 1.1);
      assert(tot_mass > 0.9);
    }
#endif
#if 0  
  for (uint32_t i = 0; i < N; i++)
    {
      Coordinates xyz;
      int32_t ixyz[NUMDIMS];
      for (int j = 0; j < NUMDIMS; j++)
	{
	  xyz[j] = (particles[i].coords[j] / spatialLen * szGrid);
	  ixyz[j] = (int32_t)round(xyz[j] - 0.5);
	  if (ixyz[j] < 0)
	    ixyz[j] = szGrid-1;
	  else if (ixyz[j] >= szGrid)
	    ixyz[j] = 0;
	}

      uint32_t idx = ixyz[0] + ixyz[1] * szGrid + ixyz[2] * szGrid * szGrid;
      densityGrid[idx] += particles[i].mass;      
    }

#endif

  for (uint32_t i = 0; i < N; i++)
    {
      CICType x, y, z;
      int32_t ix, iy, iz;
      int32_t ix2, iy2, iz2;

      x = particles[i].coords[0] / spatialLen * szGrid + 0.5;
      y = particles[i].coords[1] / spatialLen * szGrid + 0.5;
      z = particles[i].coords[2] / spatialLen * szGrid + 0.5;

      if (x < 0)
	x += szGrid;
      if (y < 0)
	y += szGrid;
      if (z < 0)
	z += szGrid;

      ix = ((int32_t)floor(x));
      iy = ((int32_t)floor(y));
      iz = ((int32_t)floor(z));
      
      ix2 = (ix + 1) % szGrid;
      iy2 = (iy + 1) % szGrid;
      iz2 = (iz + 1) % szGrid;

      CICType alpha_x = x - ix;
      CICType alpha_y = y - iy;
      CICType alpha_z = z - iz;

      ix %= szGrid;
      iy %= szGrid;
      iz %= szGrid;
      
      assert(alpha_x >= 0);
      assert(alpha_y >= 0);
      assert(alpha_z >= 0);

      CICType beta_x = 1 - alpha_x;
      CICType beta_y = 1 - alpha_y;
      CICType beta_z = 1 - alpha_z;

      assert(beta_x >= 0);
      assert(beta_y >= 0);
      assert(beta_z >= 0);

      CICType mass = particles[i].mass;
      uint32_t idx;

      // 000
      idx = ix + (iy + iz * szGrid) * szGrid;
      densityGrid[idx] +=
	mass * beta_x * beta_y * beta_z;
      
      // 100
      idx = ix2 + (iy + iz * szGrid) * szGrid;
      densityGrid[idx] +=
	mass * alpha_x * beta_y * beta_z;
      
      // 010
      idx = ix + (iy2 + iz * szGrid) * szGrid;
      densityGrid[idx] +=
	mass * beta_x * alpha_y * beta_z;
      
      // 110
      idx = ix2 + (iy2 + iz * szGrid) * szGrid;
      densityGrid[idx] +=
	mass * alpha_x * alpha_y * beta_z;

      // 001
      idx = ix + (iy + iz2 * szGrid) * szGrid;
      densityGrid[idx] +=
	mass * beta_x * beta_y * alpha_z;

      // 101
      idx = ix2 + (iy + iz2 * szGrid) * szGrid;
      densityGrid[idx] +=
	mass * alpha_x * beta_y * alpha_z;

      // 011
      idx = ix + (iy2 + iz2 * szGrid) * szGrid;
      densityGrid[idx] +=
	mass * beta_x * alpha_y * alpha_z;

      // 111
      idx = ix2 + (iy2 + iz2 * szGrid) * szGrid;
      densityGrid[idx] +=
	mass * alpha_x * alpha_y * alpha_z;
    }
}
  
void CICFilter::getDensityField(CICType*& field, uint32_t& res)
{
  field = densityGrid;
  res = totalSize;
}
