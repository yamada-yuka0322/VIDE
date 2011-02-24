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
};

bool loadParticleInfo(ParticleInfo& info,
		      const std::string& particles, 
		      const std::string& extra_info);

#endif
