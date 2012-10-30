#include <cmath>

namespace CosmoTool {

template<typename ValType, int Ndims>
SPHSmooth<ValType,Ndims>::SPHSmooth(SPHTree *tree, uint32_t Nsph)
{
  this->Nsph = Nsph;
  this->tree = tree;
  this->currentNgb = 0;

  this->maxNgb = Nsph;
  ngb = new SPHCell *[maxNgb];  
  distances = new CoordType[maxNgb];
}

template<typename ValType, int Ndims>
SPHSmooth<ValType,Ndims>::~SPHSmooth()
{
  delete[] ngb;
  delete[] distances;
}

template<typename ValType, int Ndims>
ComputePrecision SPHSmooth<ValType,Ndims>::computeWValue(const typename SPHTree::coords& c,
							 SPHCell& cell,
							 CoordType d,
							 computeParticleValue fun)
{
  CoordType weight;

  d /= smoothRadius;
  weight = getKernel(d);

  if (cell.val.weight != 0)
    return weight * fun(cell.val.pValue) / cell.val.weight;
  else
    return 0;
}

template<typename ValType, int Ndims>
void
SPHSmooth<ValType,Ndims>::fetchNeighbours(const typename SPHTree::coords& c, uint32_t newNngb)
{
  ComputePrecision d2, max_dist = 0;
  uint32_t requested = newNngb;

  if (requested > maxNgb)
    {
      delete[] ngb;
      delete[] distances;

      maxNgb = requested;
      ngb = new SPHCell *[maxNgb];
      distances = new CoordType[maxNgb];
    }

  memcpy(currentCenter, c, sizeof(c));
  tree->getNearestNeighbours(c, requested, ngb, distances);

  currentNgb = 0;
  for (uint32_t i = 0; i < requested && ngb[i] != 0; i++,currentNgb++)
    {
      distances[i] = sqrt(distances[i]);
      d2 = distances[i];
      if (d2 > max_dist)
        max_dist = d2;
    }

  smoothRadius = max_dist / 2;  
}

template<typename ValType, int Ndims>
void
SPHSmooth<ValType,Ndims>::fetchNeighbours(const typename SPHTree::coords& c)
{
  ComputePrecision d2, max_dist = 0;
  uint32_t requested = Nsph;

  memcpy(currentCenter, c, sizeof(c));
  tree->getNearestNeighbours(c, requested, ngb, distances);

  currentNgb = 0;
  for (uint32_t i = 0; i < requested && ngb[i] != 0; i++,currentNgb++)
    {
      distances[i] = sqrt(distances[i]);
      d2 = distances[i];
      if (d2 > max_dist)
	max_dist = d2;
    }

  smoothRadius = max_dist / 2;  
}					

template<typename ValType, int Ndims>
void
SPHSmooth<ValType,Ndims>::fetchNeighboursOnVolume(const typename SPHTree::coords& c,
						  ComputePrecision radius)
{
  uint32_t numPart;
  ComputePrecision d2, max_dist = 0;

  memcpy(currentCenter, c, sizeof(c));

  currentNgb = tree->getIntersection(c, radius, ngb, distances,
				  maxNgb);

  for (uint32_t i = 0; i < currentNgb; i++)
    {
      distances[i] = sqrt(distances[i]);
      d2 = distances[i];
      if (d2 > max_dist)
	max_dist = d2;
    }
  smoothRadius = max_dist / 2;  
}				

template<typename ValType, int Ndims>
ComputePrecision 
SPHSmooth<ValType,Ndims>::computeSmoothedValue(const typename SPHTree::coords& c,
					       computeParticleValue fun)
{
  ComputePrecision outputValue = 0;
  ComputePrecision max_dist = 0;
  ComputePrecision r3 = smoothRadius * smoothRadius * smoothRadius;

  for (uint32_t i = 0; i < currentNgb; i++)
    {
      outputValue += computeWValue(c, *ngb[i], distances[i], fun);
    }

  return outputValue / r3;
}

template<typename ValType>
ComputePrecision interpolateOne(const ValType& t)
{
  return 1.0;
}

// WARNING ! Cell's weight must be 1 !!!
template<typename ValType, int Ndims>
ComputePrecision SPHSmooth<ValType,Ndims>::computeInterpolatedValue(const typename SPHTree::coords& c,
								    computeParticleValue fun)
{
  ComputePrecision outputValue = 0;
  ComputePrecision max_dist = 0;
  ComputePrecision weight = 0;

  for (uint32_t i = 0; i < currentNgb; i++)
    {
      outputValue += computeWValue(c, *ngb[i], distances[i], fun);
      weight += computeWValue(c, *ngb[i], distances[i], interpolateOne);
    }

  return (outputValue == 0) ? 0 : (outputValue / weight);  
}

template<typename ValType, int Ndims>
void SPHSmooth<ValType,Ndims>::runForEachNeighbour(runParticleValue fun)
{
  for (uint32_t i = 0; i < currentNgb; i++)
    {
      fun(ngb[i]);
    }
}


template<typename ValType, int Ndims>
void SPHSmooth<ValType,Ndims>::addGridSite(const typename SPHTree::coords& c)
{
  ComputePrecision outputValue = 0;
  ComputePrecision max_dist = 0;
  
  ComputePrecision r3 = smoothRadius * smoothRadius * smoothRadius;

  for (uint32_t i = 0; i < currentNgb; i++)
    {
      ComputePrecision d = distances[i];
      SPHCell& cell = *ngb[i];
      cell.val.weight += getKernel(d/smoothRadius) / r3;

    }
}

template<typename ValType, int Ndims>
ComputePrecision 
SPHSmooth<ValType,Ndims>::getKernel(ComputePrecision x) const
{
  // WARNING !!! This is an unnormalized version of the kernel.
  if (x < 1)
    return 1 - 1.5 * x * x + 0.75 * x * x * x;
  else if (x < 2)
    {
      ComputePrecision d = 2 - x;
      return 0.25 * d * d * d;
    }
  else
    return 0;
}

template<typename ValType, int Ndims>
bool SPHSmooth<ValType,Ndims>::hasNeighbours() const
{
  return (currentNgb != 0);
}

template<class ValType1, class ValType2, int Ndims>
bool operator<(const SPHSmooth<ValType1, Ndims>& s1, const SPHSmooth<ValType2, Ndims>& s2)
{
  return (s1.getSmoothingLen() < s2.getSmoothingLen());
}

};
