#ifndef __CICFILTER_HPP
#define __CICFILTER_HPP

#include "CosmoTool/config.hpp"
#include <inttypes.h>

using namespace CosmoTool;

typedef float CICType;

  typedef struct
  {
    float mass;
    Coordinates coords;
  } CICParticles;

  class CICFilter
  {
  public:
    CICFilter(uint32_t resolution, double spatialLen);
    ~CICFilter();

    void resetMesh();
    void putParticles(CICParticles *particles, uint32_t N);    

    void getDensityField(CICType*& field, uint32_t& res);

  protected:
    CICType *densityGrid;
    double spatialLen;
    uint32_t totalSize;
    uint32_t szGrid;
  };

#endif
