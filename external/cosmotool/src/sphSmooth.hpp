#ifndef __COSMOTOOL_SPH_SMOOTH_HPP
#define __COSMOTOOL_SPH_SMOOTH_HPP

#include "config.hpp"
#include "mykdtree.hpp"

namespace CosmoTool
{
  template <typename ValType, int Ndims = NUMDIMS>
  class SPHSmooth
  {
  public:
    typedef struct
    {
      ComputePrecision weight;
      ValType pValue;
    } FullType;

    typedef KDTree<Ndims,FullType> SPHTree;
    typedef KDTreeNode<Ndims,FullType> SPHNode;
    typedef KDCell<Ndims,FullType>  SPHCell;
    typedef typename KDTree<Ndims,FullType>::CoordType CoordType;

    typedef ComputePrecision (*computeParticleValue)(const ValType& t);
    typedef void (*runParticleValue)(ValType& t);

  public:
    SPHSmooth(SPHTree *tree, uint32_t Nsph);
    virtual ~SPHSmooth();   

    void fetchNeighbours(const typename SPHTree::coords& c);
    void fetchNeighbours(const typename SPHTree::coords& c, uint32_t newNsph);
    void fetchNeighboursOnVolume(const typename SPHTree::coords& c, ComputePrecision radius);
    const typename SPHTree::coords& getCurrentCenter() const
    {
      return currentCenter;
    }

    ComputePrecision computeSmoothedValue(const typename SPHTree::coords& c,
					  computeParticleValue fun);    
    ComputePrecision computeInterpolatedValue(const typename SPHTree::coords& c,
					      computeParticleValue fun);    
    ComputePrecision getMaxDistance(const typename SPHTree::coords& c, 
				    SPHNode *node) const;

    ComputePrecision getSmoothingLen() const
    {
      return smoothRadius;
    }

    // TO USE WITH EXTREME CARE !
    void setSmoothingLen(ComputePrecision len)
    {
      smoothRadius = len;
    }
    // END

    void runForEachNeighbour(runParticleValue fun);
    void addGridSite(const typename SPHTree::coords& c);

    bool hasNeighbours() const;

    virtual ComputePrecision getKernel(ComputePrecision d) const;    

    SPHTree *getTree() { return tree; }

    void changeNgb(uint32_t newMax) {
      delete[] ngb;
      delete[] distances;
      ngb = new SPHCell *[newMax];
      distances = new CoordType[newMax];
      maxNgb = newMax;
    }

    uint32_t getCurrent() const { return currentNgb; }

    uint32_t getNgb() const { return maxNgb; }
 
  protected:
    SPHCell **ngb;
    CoordType *distances;
    uint32_t Nsph;
    uint32_t deltaNsph;
    uint32_t maxNgb;
    uint32_t currentNgb;
    SPHTree *tree;
    ComputePrecision smoothRadius;
    typename SPHTree::coords currentCenter;
    
    ComputePrecision computeWValue(const typename SPHTree::coords & c,
				   SPHCell& cell,
				   CoordType d,
				   computeParticleValue fun);
    void runUnrollNode(SPHNode *node,
		       runParticleValue fun);
  };

  template<class ValType1, class ValType2, int Ndims>
  bool operator<(const SPHSmooth<ValType1, Ndims>& s1, const SPHSmooth<ValType2, Ndims>& s2);

};

#include "sphSmooth.tcc"

#endif
