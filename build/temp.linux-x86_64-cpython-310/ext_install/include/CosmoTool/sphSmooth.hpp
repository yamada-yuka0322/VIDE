/*+
This is CosmoTool (./src/sphSmooth.hpp) -- Copyright (C) Guilhem Lavaux (2007-2014)

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

#ifndef __COSMOTOOL_SPH_SMOOTH_HPP
#define __COSMOTOOL_SPH_SMOOTH_HPP
#include <boost/shared_ptr.hpp>
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
    typedef SPHCell *P_SPHCell;
    struct SPHState
    {
      boost::shared_ptr<P_SPHCell[]> ngb;
      boost::shared_ptr<CoordType[]> distances;
      typename SPHTree::coords currentCenter;
      int currentNgb;
      ComputePrecision smoothRadius;
    };
    
  
    SPHSmooth(SPHTree *tree, uint32_t Nsph);
    virtual ~SPHSmooth();   

    void fetchNeighbours(const typename SPHTree::coords& c, SPHState *state = 0);
    
    void fetchNeighbours(const typename SPHTree::coords& c, uint32_t newNsph);
    void fetchNeighboursOnVolume(const typename SPHTree::coords& c, ComputePrecision radius);
    const typename SPHTree::coords& getCurrentCenter() const
    {
      return internal.currentCenter;
    }

    template<typename FuncT>
    ComputePrecision computeSmoothedValue(const typename SPHTree::coords& c,
					  FuncT fun, SPHState *state = 0);    
	
	template<typename FuncT>
    ComputePrecision computeInterpolatedValue(const typename SPHTree::coords& c,
					      FuncT fun, SPHState *state = 0);
					          
    ComputePrecision getMaxDistance(const typename SPHTree::coords& c, 
				    SPHNode *node) const;

    ComputePrecision getSmoothingLen() const
    {
      return internal.smoothRadius;
    }

    // TO USE WITH EXTREME CARE !
    void setSmoothingLen(ComputePrecision len)
    {
      internal.smoothRadius = len;
    }
    // END

    template<typename FuncT>
    void runForEachNeighbour(FuncT fun, SPHState *state = 0);
    void addGridSite(const typename SPHTree::coords& c);

    bool hasNeighbours() const;

    virtual ComputePrecision getKernel(ComputePrecision d) const;    

    SPHTree *getTree() { return tree; }

    void changeNgb(uint32_t newMax) {
      internal.ngb = boost::shared_ptr<SPHCell *>(new SPHCell *[newMax]);
      internal.distances = boost::shared_ptr<CoordType>(new CoordType[newMax]);
      maxNgb = newMax;
    }

    uint32_t getCurrent() const { return internal.currentNgb; }

    uint32_t getNgb() const { return maxNgb; }
 
  protected:
    SPHState internal;
    uint32_t Nsph;
    uint32_t deltaNsph;
    uint32_t maxNgb;
    SPHTree *tree;
    
    template<typename FuncT>
    ComputePrecision computeWValue(const typename SPHTree::coords & c,
				   SPHCell& cell,
				   CoordType d,
				   FuncT fun, SPHState *state);
    
    template<typename FuncT>
    void runUnrollNode(SPHNode *node,
		       FuncT fun);
  };

  template<class ValType1, class ValType2, int Ndims>
  bool operator<(const SPHSmooth<ValType1, Ndims>& s1, const SPHSmooth<ValType2, Ndims>& s2);

};

#include "sphSmooth.tcc"

#endif
