/*+
This is CosmoTool (./src/octTree.hpp) -- Copyright (C) Guilhem Lavaux (2007-2014)

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

#ifndef __COSMOTOOL_AMR_HPP
#define __COSMOTOOL_AMR_HPP

#include "config.hpp"

namespace CosmoTool
{

  typedef uint32_t octPtr;
  typedef uint16_t octCoordType;
  static const uint16_t octCoordTypeNorm = 0xffff;
  static const uint16_t octCoordCenter = 0x8000;

  // This is also the root cell, but this one
  // is never referenced (really ??).
  static const octPtr invalidOctCell = 0xffffffff;
  static const octPtr emptyOctCell = 0;
  static const octPtr octParticleMarker = 0x80000000;
  static const octPtr octParticleMask = 0x7fffffff;

  typedef octCoordType OctCoords[3];

  template<class T = void>
  struct OctCell
  {
    octPtr numberLeaves;
    octPtr children[8];
    T data;
  };

  template<typename T>
  struct OctTree_defaultUpdater {
    void operator()(T& d) { }
  };

  template<typename T_dataUpdater = OctTree_defaultUpdater<void>, class T = void>
  class OctTree
  {
  public:
    OctTree(const FCoordinates *particles, octPtr numParticles,
	    uint32_t maxTreeDepth,  uint32_t maxAbsoluteDepth,
	    uint32_t threshold = 1);
    ~OctTree();

    void buildTree(uint32_t maxAbsoluteDepth);
    void insertParticle(octPtr node, 
			const OctCoords& icoord,
			octCoordType halfNodeLength,
			octPtr particleId,
			uint32_t maxAbsoluteDepth);


    octPtr getNumberLeaves() const {
      return cells[0].numberLeaves;
    }

    static bool unconditioned(const FCoordinates&, octPtr, float, bool)
    {
      return true;
    }

    template<typename FunT>
    void walkTree(FunT f)
    {
      walkTree(f, unconditioned);
    }

    template<typename FunT, typename CondT>
    void walkTree(FunT f, CondT condition)
    {
      OctCoords rootCenter = { octCoordCenter, octCoordCenter, octCoordCenter };

      walkTreeElements(f, condition, 0, rootCenter, octCoordCenter);
    }

    

  protected:
    T_dataUpdater updater;
    const FCoordinates *particles;
    octPtr numParticles;
    OctCell<T> *cells;
    float Lbox;
    octPtr lastNode;
    octPtr numCells;  
    float lenNorm;
    float xMin[3];


    static bool unconditioned()
    {
      return true;
    }

    template<typename FunT, typename CondT>
    void walkTreeElements(FunT f, CondT condition,
			  octPtr node, 
			  const OctCoords& icoord,
			  octCoordType halfNodeLength)
    {
      OctCoords newCoord;
      FCoordinates center, realCenter;

      for (int j = 0; j < 3; j++)
        {
          center[j] = icoord[j]/(2.*octCoordCenter);
          realCenter[j] = xMin[j] + center[j]*lenNorm;
        }

      f(realCenter, cells[node].numberLeaves, lenNorm*halfNodeLength/(float)octCoordCenter,
        cells[node].children[0] == invalidOctCell, // True if this is a meta-node
        false);
     
      if (!condition(realCenter, cells[node].numberLeaves,
		     lenNorm*halfNodeLength/(float)octCoordCenter,
		     cells[node].children[0] == invalidOctCell))
        return;
       
      for (int i = 0; i < 8; i++)
        {
          octPtr newNode = cells[node].children[i];
          int ipos[3] = { (i&1), (i&2)>>1, (i&4)>>2 };
          
          if (newNode == emptyOctCell || newNode == invalidOctCell)
            continue;
          	  
          for (int j = 0; j < 3; j++)
            newCoord[j] = icoord[j]+(2*ipos[j]-1)*halfNodeLength/2;
          
          if (newNode & octParticleMarker)
            {
              for (int j = 0; j < 3; j++)
                {
                  center[j] = newCoord[j]/(2.*octCoordCenter);
                  realCenter[j] = xMin[j] + lenNorm*center[j];
                }

              f(realCenter, 
                1, lenNorm*halfNodeLength/(2.*octCoordCenter),
                false, true);
              continue;
            }

          walkTreeElements(f, condition, cells[node].children[i], newCoord, halfNodeLength/2);
        }
      
    }

  };
  
};


#include "octTree.tcc"

#endif
