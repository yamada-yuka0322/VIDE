/*+
This is CosmoTool (./src/octTree.cpp) -- Copyright (C) Guilhem Lavaux (2007-2014)

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

#include <iostream>
#include <cmath>
#include <cassert>
#include "config.hpp"
#include "octTree.hpp"

namespace CosmoTool {

using namespace std;

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

template<typename Updater, typename T>
OctTree<Updater,T>::OctTree(const FCoordinates *particles, octPtr numParticles, 
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

  cells = new OctCell<T>[numCells];
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

template<typename Updater, typename T>
OctTree<Updater,T>::~OctTree()
{
  delete cells;
}

template<typename Updater, typename T>
void OctTree<Updater,T>::buildTree(uint32_t maxAbsoluteDepth)
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


template<typename Updater, typename T>
void OctTree<Updater,T>::insertParticle(octPtr node, 
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


};
