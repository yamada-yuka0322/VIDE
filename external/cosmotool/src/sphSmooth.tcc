#include <cmath>
#include "algo.hpp"

namespace CosmoTool {

template<typename ValType, int Ndims>
SPHSmooth<ValType,Ndims>::SPHSmooth(SPHTree *tree, uint32_t Nsph)
{
  this->Nsph = Nsph;
  this->tree = tree;
  internal.currentNgb = 0;

  this->maxNgb = Nsph;
  internal.ngb = boost::shared_ptr<SPHCell *[]>(new SPHCell *[maxNgb]);  
  internal.distances = boost::shared_ptr<CoordType[]>(new CoordType[maxNgb]);
}

template<typename ValType, int Ndims>
SPHSmooth<ValType,Ndims>::~SPHSmooth()
{
}

template<typename ValType, int Ndims>
template<typename FuncT>
ComputePrecision SPHSmooth<ValType,Ndims>::computeWValue(const typename SPHTree::coords& c,
							 SPHCell& cell,
							 CoordType d,
							 FuncT fun, SPHState *state)
{
  CoordType weight;

  d /= state->smoothRadius;
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
      maxNgb = requested;
      internal.ngb = boost::shared_ptr<P_SPHCell[]>(new P_SPHCell[maxNgb]);
      internal.distances = boost::shared_ptr<CoordType[]>(new CoordType[maxNgb]);
    }

  memcpy(internal.currentCenter, c, sizeof(c));
  tree->getNearestNeighbours(c, requested, (SPHCell **)internal.ngb.get(), (CoordType*)internal.distances.get());

  internal.currentNgb = 0;
  for (uint32_t i = 0; i < requested && (internal.ngb)[i] != 0; i++,internal.currentNgb++)
    {
      internal.distances[i] = sqrt(internal.distances[i]);
      d2 = internal.distances[i];
      if (d2 > max_dist)
        max_dist = d2;
    }

  internal.smoothRadius = max_dist / 2;  
}

template<typename ValType, int Ndims>
void SPHSmooth<ValType,Ndims>::fetchNeighbours(const typename SPHTree::coords& c, SPHState *state)
{
  ComputePrecision d2, max_dist = 0;
  uint32_t requested = Nsph;

  if (state != 0) {
    state->distances = boost::shared_ptr<CoordType[]>(new CoordType[Nsph]);
    state->ngb = boost::shared_ptr<SPHCell *[]>(new SPHCell *[Nsph]);
  } else
    state = &internal;
    
  memcpy(state->currentCenter, c, sizeof(c));
  
  tree->getNearestNeighbours(c, requested, state->ngb.get(), state->distances.get());

  state->currentNgb = 0;
  for (uint32_t i = 0; i < requested && state->ngb[i] != 0; i++,state->currentNgb++)
    {
      d2 = internal.distances[i] = sqrt(internal.distances[i]);
      if (d2 > max_dist)
        max_dist = d2;
    }

  state->smoothRadius = max_dist / 2;  
}


template<typename ValType, int Ndims>
void
SPHSmooth<ValType,Ndims>::fetchNeighboursOnVolume(const typename SPHTree::coords& c,
                                                  ComputePrecision radius)
{
  uint32_t numPart;
  ComputePrecision d2, max_dist = 0;

  memcpy(internal.currentCenter, c, sizeof(c));

  internal.currentNgb = tree->getIntersection(c, radius, internal.ngb, internal.distances,
				  maxNgb);

  for (uint32_t i = 0; i < internal.currentNgb; i++)
    {
      d2 = internal.distances[i] = sqrt(internal.distances[i]);
      if (d2 > max_dist)
        max_dist = d2;
    }
  internal.smoothRadius = max_dist / 2;  
}				

template<typename ValType, int Ndims>
template<typename FuncT>
ComputePrecision 
SPHSmooth<ValType,Ndims>::computeSmoothedValue(const typename SPHTree::coords& c,
                                               FuncT fun, SPHState *state)
{
  if (state == 0)
    state = &internal;
    
  ComputePrecision outputValue = 0;
  ComputePrecision max_dist = 0;
  ComputePrecision r3 = cube(state->smoothRadius);

  for (uint32_t i = 0; i < state->currentNgb; i++)
    {
      outputValue += computeWValue(c, *state->ngb[i], state->distances[i], fun, state);
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
template<typename FuncT>
ComputePrecision SPHSmooth<ValType,Ndims>::computeInterpolatedValue(const typename SPHTree::coords& c,
								    FuncT fun, SPHState *state)
{
  if (state == 0)
    state = &internal;
    
  ComputePrecision outputValue = 0;
  ComputePrecision max_dist = 0;
  ComputePrecision weight = 0;

  for (uint32_t i = 0; i < state->currentNgb; i++)
    {
      outputValue += computeWValue(c, *state->ngb[i], state->distances[i], fun);
      weight += computeWValue(c, *state->ngb[i], state->distances[i], interpolateOne);
    }

  return (outputValue == 0) ? 0 : (outputValue / weight);  
}

template<typename ValType, int Ndims>
template<typename FuncT>
void SPHSmooth<ValType,Ndims>::runForEachNeighbour(FuncT fun, SPHState *state)
{
  if (state == 0)
    state = &internal;
    
  for (uint32_t i = 0; i < state->currentNgb; i++)
    {
      fun(state->ngb[i]);
    }
}


template<typename ValType, int Ndims>
void SPHSmooth<ValType,Ndims>::addGridSite(const typename SPHTree::coords& c)
{
  ComputePrecision outputValue = 0;
  ComputePrecision max_dist = 0;
  
  ComputePrecision r3 = cube(internal.smoothRadius);

  for (uint32_t i = 0; i < internal.currentNgb; i++)
    {
      ComputePrecision d = internal.distances[i];
      SPHCell& cell = *(internal.ngb[i]);
      cell.val.weight += getKernel(d/internal.smoothRadius) / r3;
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
  return (internal.currentNgb != 0);
}

template<class ValType1, class ValType2, int Ndims>
bool operator<(const SPHSmooth<ValType1, Ndims>& s1, const SPHSmooth<ValType2, Ndims>& s2)
{
  return (s1.getSmoothingLen() < s2.getSmoothingLen());
}

};
