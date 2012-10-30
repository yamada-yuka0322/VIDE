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

  struct OctCell
  {
    octPtr numberLeaves;
    octPtr children[8];
  };

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
    const FCoordinates *particles;
    octPtr numParticles;
    OctCell *cells;
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


#endif
