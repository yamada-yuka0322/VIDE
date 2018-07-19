/*+
This is CosmoTool (./src/cic.cpp) -- Copyright (C) Guilhem Lavaux (2007-2014)

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
#include "openmp.hpp"
#include <assert.h>
#include <math.h>
#include <inttypes.h>
#include "cic.hpp"

using namespace CosmoTool;

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
  int threadUsed = smp_get_max_threads();
  double *threadedDensity[threadUsed];
  bool threadActivated[threadUsed];
  uint32_t tUsedMin[threadUsed], tUsedMax[threadUsed];

  for (int t = 0; t < threadUsed; t++)
  {
    threadedDensity[t] = new double[totalSize];
    std::fill(threadedDensity[t], threadedDensity[t]+totalSize, 0);
  }

  std::fill(threadActivated, threadActivated+threadUsed, false);
  std::fill(tUsedMin, tUsedMin+threadUsed, totalSize);
  std::fill(tUsedMax, tUsedMax+threadUsed, 0);

#pragma omp parallel
  {
    int thisThread = smp_get_thread_id();
    double *dg = threadedDensity[thisThread];

    threadActivated[thisThread] = true;
  
#pragma omp for schedule(static)
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
      dg[idx] +=
	mass * beta_x * beta_y * beta_z;
      
      // 100
      idx = ix2 + (iy + iz * szGrid) * szGrid;
      dg[idx] +=
	mass * alpha_x * beta_y * beta_z;
      
      // 010
      idx = ix + (iy2 + iz * szGrid) * szGrid;
      dg[idx] +=
	mass * beta_x * alpha_y * beta_z;
      
      // 110
      idx = ix2 + (iy2 + iz * szGrid) * szGrid;
      dg[idx] +=
	mass * alpha_x * alpha_y * beta_z;

      // 001
      idx = ix + (iy + iz2 * szGrid) * szGrid;
      dg[idx] +=
	mass * beta_x * beta_y * alpha_z;

      // 101
      idx = ix2 + (iy + iz2 * szGrid) * szGrid;
      dg[idx] +=
	mass * alpha_x * beta_y * alpha_z;

      // 011
      idx = ix + (iy2 + iz2 * szGrid) * szGrid;
      dg[idx] +=
	mass * beta_x * alpha_y * alpha_z;

      // 111
      idx = ix2 + (iy2 + iz2 * szGrid) * szGrid;
      dg[idx] +=
	mass * alpha_x * alpha_y * alpha_z;

      tUsedMin[thisThread] = std::min(tUsedMin[thisThread], idx);
      tUsedMax[thisThread] = std::max(tUsedMax[thisThread], idx);
    }
  }

  for (int t = 0; t < threadUsed; t++)
  {
    if (!threadActivated[t])
      continue;
#pragma omp parallel for schedule(static)
    for (long p = tUsedMin[t]; p < tUsedMax[t]; p++)
      densityGrid[p] += threadedDensity[t][p];
  }

  for (int t = 0; t < threadUsed; t++)
  {
    delete[] threadedDensity[t];
  }
}
  
void CICFilter::getDensityField(CICType*& field, uint32_t& res)
{
  field = densityGrid;
  res = totalSize;
}

