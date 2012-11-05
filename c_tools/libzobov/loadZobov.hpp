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
