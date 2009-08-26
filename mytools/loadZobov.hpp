#ifndef __LOAD_ZOBOV_HPP
#define __LOAD_ZOBOV_HPP

#include <vector>

struct ZobovZone
{
  std::vector<int> pId;
};

struct ZobovVoid
{
  std::vector<ZobovZone *> zId;
  float proba;
  int numParticles, coreParticle;
  float volume;
};

struct ZobovRep
{
  std::vector<ZobovZone> allZones;
  std::vector<ZobovVoid> allVoids;
};

bool loadZobov(const char *descName,
	       const char *adjName, const char *volName, ZobovRep& z);

#endif
