#include <iostream>
#include <cmath>
#include <cassert>
#include "config.hpp"
#include "octTree.hpp"

using namespace std;
using namespace CosmoTool;

//#define VERBOSE

static uint32_t mypow(uint32_t i, uint32_t p)
{
  if (p == 0)
    return 1;
  else if (p == 1)
    return i;

  uint32_t k = p/2;
  uint32_t j = mypow(i, k);
  if (2*k==p)
    return j*j;
  else
    return j*j*i;  
}

OctTree::OctTree(const FCoordinates *particles, octPtr numParticles, 
		 uint32_t maxMeanTreeDepth, uint32_t maxAbsoluteDepth,
		 uint32_t threshold)
{
  cout << "MeanTree=" << maxMeanTreeDepth << endl;
  numCells = mypow(8, maxMeanTreeDepth);
  assert(numCells < invalidOctCell);
  //#ifdef VERBOSE
  cerr << "Allocating " << numCells << " octtree cells" << endl;
  //#endif
  
  for (int j = 0; j < 3; j++)
    xMin[j] = particles[0][j];

  for (octPtr i = 1; i < numParticles; i++)
    {
      for (int j = 0; j < 3; j++)
	{
	  if (particles[i][j] < xMin[j])
	    xMin[j] = particles[i][j];      
	}
    }
  
  lenNorm = 0;
  for (octPtr i = 0; i < numParticles; i++)
    {
      for (int j = 0; j < 3; j++)
	{
	  float delta = particles[i][j]-xMin[j];
	  if (delta > lenNorm)
	    lenNorm = delta;
	}
    }
  cout << xMin[0] << " " << xMin[1] << " " << xMin[2] << " lNorm=" << lenNorm << endl; 

  cells = new OctCell[numCells];
  Lbox = (float)(octCoordTypeNorm+1);

  cells[0].numberLeaves = 0;
  for (int i = 0; i < 8; i++)
    cells[0].children[i] = emptyOctCell;

  lastNode = 1;
  this->particles = particles;
  this->numParticles = numParticles;
  buildTree(maxAbsoluteDepth);
  //#ifdef VERBOSE
  cerr << "Used " << lastNode << " cells" << endl;
  //#endif
}

OctTree::~OctTree()
{
  delete cells;
}

void OctTree::buildTree(uint32_t maxAbsoluteDepth)
{
  for (octPtr i = 0; i < numParticles; i++)
    {
      OctCoords rootCenter = { octCoordCenter, octCoordCenter, octCoordCenter };
      insertParticle(0, // root node
		     rootCenter,
		     octCoordCenter,
		     i,
		     maxAbsoluteDepth);
    }
}


void OctTree::insertParticle(octPtr node, 
			     const OctCoords& icoord,
			     octCoordType halfNodeLength,
			     octPtr particleId,
			     uint32_t maxAbsoluteDepth)
{
  
#ifdef VERBOSE
  cout << "Entering " << node << " (" << icoord[0] << "," << icoord[1] << "," << icoord[2] << ")" << endl;
#endif
  int octPos = 0;
  int ipos[3] = { 0,0,0};
  octPtr newNode;
  OctCoords newCoord;
  
  cells[node].numberLeaves++;
  if (maxAbsoluteDepth == 0)
    {
      // All children must be invalid.
      for (int i = 0 ; i < 8; i++)
	cells[node].children[i] = invalidOctCell;
      
      return;
    }
  
  for (int j = 0; j < 3; j++)
    {
      float treePos = (particles[particleId][j]-xMin[j])*Lbox/lenNorm;
      if ((octPtr)(treePos) > icoord[j])
	{
	  octPos |= (1 << j);
	  ipos[j] = 1;
	}
    }
  
  if (cells[node].children[octPos] == emptyOctCell)
    {
      // Put the particle there.
      cells[node].children[octPos] = particleId | octParticleMarker;
      return;
    }
  
  // If it is a node, explores it.
  if (!(cells[node].children[octPos] & octParticleMarker))
    {
      assert(halfNodeLength >= 2);
      // Compute coordinates
      for (int j = 0; j < 3; j++)
	newCoord[j] = icoord[j]+(2*ipos[j]-1)*halfNodeLength/2;
      insertParticle(cells[node].children[octPos], newCoord, halfNodeLength/2, 
		     particleId, maxAbsoluteDepth-1);
      return;
    }
  
  // We have a particle there. 
  // Make a new node and insert the old particle into this node.
  // Insert the new particle into the node also
  // Finally put the node in place
  
  newNode = lastNode++;
  assert(lastNode != numCells);
  
  for (int j = 0; j < 8; j++)
    cells[newNode].children[j] = emptyOctCell;
  cells[newNode].numberLeaves = 0;
  
  // Compute coordinates
  for (int j = 0; j < 3; j++)
    newCoord[j] = icoord[j]+(2*ipos[j]-1)*halfNodeLength/2;
  
  octPtr oldPartId = cells[node].children[octPos] & octParticleMask;
  
  insertParticle(newNode, newCoord, halfNodeLength/2,
		 oldPartId, maxAbsoluteDepth-1);
  insertParticle(newNode, newCoord, halfNodeLength/2,
		 particleId, maxAbsoluteDepth-1);  
  cells[node].children[octPos] = newNode;
}

