/*+
This is CosmoTool (./src/kdtree_leaf.hpp) -- Copyright (C) Guilhem Lavaux (2007-2014)

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

#ifndef __LEAF_KDTREE_HPP
#define __LEAF_KDTREE_HPP

#include <cmath>
#include "config.hpp"
#include "bqueue.hpp"

namespace CosmoTool {

  template<int N, typename CType = ComputePrecision> 
  struct KDLeafDef
  {
    typedef CType CoordType;
    typedef float KDLeafCoordinates[N];
  };
  
  template<int N, typename ValType, typename CType = ComputePrecision>
  struct KDLeafCell
  {
    bool active;
    ValType val;
    typename KDLeafDef<N,CType>::KDLeafCoordinates coord;
  };

  class NotEnoughCells: public Exception
  {
  public:
    NotEnoughCells() : Exception() {}
    ~NotEnoughCells() throw () {}
  };

  template<int N, typename ValType, typename CType = ComputePrecision>
  struct KDLeafTreeNode
  {
    bool leaf;
    union {
      KDLeafCell<N,ValType,CType> *value;
      KDLeafTreeNode<N,ValType,CType> *children[2];
    };
    typename KDLeafDef<N,CType>::KDLeafCoordinates minBound, maxBound;
#ifdef __KDLEAF_TREE_NUMNODES
    uint32_t numNodes;
#endif
  };

  template<int N, typename ValType, typename CType = ComputePrecision>
  class KDLeafTree
  {
  public:
    typedef typename KDLeafDef<N,CType>::CoordType CoordType;
    typedef typename KDLeafDef<N>::KDLeafCoordinates coords;
    typedef KDLeafCell<N,ValType,CType> Cell;
    typedef KDLeafTreeNode<N,ValType,CType> Node;

    KDLeafTree(Cell *cells, uint32_t Ncells);
    ~KDLeafTree();

    Node *getRoot() { return root; }

    void optimize();

    Node *getAllNodes() { return nodes; }
    uint32_t getNumNodes() const { return lastNode; }

    uint32_t countActives() const;

    CoordType computeDistance(const Cell *cell, const coords& x) const;

#ifdef __KDLEAF_TREE_NUMNODES
    uint32_t getNumberInNode(const Node *n) const { return n->numNodes; }
#else
    uint32_t getNumberInNode(const Node *n) const {
      if (n->leaf)
        return 1;
      if (n == 0) 
        return 0;
      return getNumberInNode(n->children[0])+getNumberInNode(n->children[1]);
    }
#endif

    double countInRange(CType sLo, CType sHi, Node *root1 = 0, Node *root2 = 0) const;

  protected:
    Node *nodes;
    uint32_t numNodes, numCells;
    uint32_t lastNode;

    Node *root;
    Cell **sortingHelper;

    Node *buildTree(Cell **cell0,
		    uint32_t NumCells,
		    uint32_t depth,
		    coords minBound,
		    coords maxBound);
    
    double recursiveCountInRange(Node *na, Node *nb, CType sLo, CType sHi) const;
  };

  template<int N, typename T, typename CType>
  uint32_t gatherActiveCells(KDLeafCell<N,T,CType> **cells, uint32_t numCells);

};

#include "kdtree_leaf.tcc"

#endif
