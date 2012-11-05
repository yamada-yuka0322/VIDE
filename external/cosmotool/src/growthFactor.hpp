#ifndef COSMO_GROWTH_FACTOR_HPP
#define COSMO_GROWTH_FACTOR_HPP

#include "interpolate.hpp"

namespace CosmoTool
{
  Interpolate buildLinearGrowth(double OmegaLambda, double OmegaMatter, double Hubble, int numPoints = 10000);

};

#endif
