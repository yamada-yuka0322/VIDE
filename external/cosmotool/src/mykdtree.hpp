/*+
This is CosmoTool (./src/mykdtree.hpp) -- Copyright (C) Guilhem Lavaux (2007-2014)

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

#ifndef __HV_KDTREE_HPP
#define __HV_KDTREE_HPP

#include <cmath>
#include "config.hpp"
#include "bqueue.hpp"

#ifndef __KD_TREE_ACTIVE_CELLS
#define __KD_TREE_ACTIVE_CELLS 1
#endif

namespace CosmoTool {

  typedef uint64_t NodeIntType;

  template<int N, typename CType = ComputePrecision> 
  struct KDDef
  {
    typedef CType CoordType;
    typedef float KDCoordinates[N];
  };
  
  template<int N, typename ValType, typename CType = ComputePrecision>
  struct KDCell
  {
#if __KD_TREE_ACTIVE_CELLS == 1
    bool active;
#endif
    ValType val;
    typename KDDef<N,CType>::KDCoordinates coord;
  };

  class NotEnoughCells: public Exception
  {
  public:
    NotEnoughCells() : Exception() {}
    ~NotEnoughCells() throw () {}
  };

  class InvalidOnDiskKDTree : public Exception
  {
  public:
    InvalidOnDiskKDTree() : Exception() {}
    ~InvalidOnDiskKDTree() throw () {}
  };

  template<int N, typename ValType, typename CType = ComputePrecision>
  struct KDTreeNode
  {
    KDCell<N,ValType,CType> *value;
    KDTreeNode<N,ValType,CType> *children[2];
    typename KDDef<N,CType>::KDCoordinates minBound, maxBound;
#ifdef __KD_TREE_NUMNODES
    NodeIntType numNodes;
#endif
  };

  template<int N, typename ValType, typename CType = ComputePrecision>
  class RecursionInfoCells
  {
  public:

    typename KDDef<N,CType>::KDCoordinates x;
    typename KDDef<N,CType>::CoordType r, r2;
    KDCell<N, ValType,CType> **cells;
    typename KDDef<N,CType>::CoordType *distances;
    uint32_t currentRank;
    uint32_t numCells;
  };
  

  template<int N, typename ValType, typename CType = ComputePrecision>
  class RecursionMultipleInfo
  {
  public:
    typename KDDef<N,CType>::KDCoordinates x;
    BoundedQueue< KDCell<N,ValType,CType> *, typename KDDef<N,CType>::CoordType> queue;
    int traversed;

    RecursionMultipleInfo(const typename KDDef<N,CType>::KDCoordinates& rx,
			  KDCell<N,ValType,CType> **cells,
			  uint32_t numCells)
      : queue(cells, numCells, INFINITY),traversed(0)
    {
      std::copy(rx, rx+N, x);
    }
  };

  template<int N, typename ValType, typename CType = ComputePrecision> 
  struct KD_default_cell_splitter
  {
    void operator()(KDCell<N,ValType,CType> **cells, NodeIntType Ncells, NodeIntType& split_index, int axis, typename KDDef<N,CType>::KDCoordinates minBound, typename KDDef<N,CType>::KDCoordinates maxBound);
  };

  template<int N, typename ValType, typename CType = ComputePrecision, typename CellSplitter = KD_default_cell_splitter<N,ValType,CType> >
  class KDTree
  {
  public:
    typedef typename KDDef<N,CType>::CoordType CoordType;
    typedef typename KDDef<N>::KDCoordinates coords;
    typedef KDCell<N,ValType,CType> Cell;
    typedef KDTreeNode<N,ValType,CType> Node;
    
    CellSplitter splitter;

    KDTree(Cell *cells, NodeIntType Ncells);
    ~KDTree();

    void setPeriodic(bool on, CoordType replicate)
     {
       periodic = on;
       std::fill(this->replicate, this->replicate+N, replicate);
     }

    void setPeriodic(bool on, const coords& replicate)
     {
       periodic = on;
       std::copy(replicate, replicate+N, this->replicate);
     }

    uint32_t getIntersection(const coords& x, CoordType r, 
			     Cell **cells,
			     uint32_t numCells);
    uint32_t getIntersection(const coords& x, CoordType r, 
			     Cell **cells,
			     CoordType *distances,
			     uint32_t numCells);
    uint32_t countCells(const coords& x, CoordType r);

    Cell *getNearestNeighbour(const coords& x);

    void getNearestNeighbours(const coords& x, uint32_t NumCells,
			      Cell **cells);
    void getNearestNeighbours(const coords& x, uint32_t NumCells,
			      Cell **cells,
			      CoordType *distances);

    Node *getRoot() { return root; }

    void optimize();

    Node *getAllNodes() { return nodes; }
    NodeIntType getNumNodes() const { return lastNode; }

    NodeIntType countActives() const;

#ifdef __KD_TREE_NUMNODES
    NodeIntType getNumberInNode(const Node *n) const { return n->numNodes; }
#else
    NodeIntType getNumberInNode(const Node *n) const {
      if (n == 0) 
        return 0;
      return 1+getNumberInNode(n->children[0])+getNumberInNode(n->children[1]);
    }
#endif

#ifdef __KD_TREE_SAVE_ON_DISK
    KDTree(std::istream& i, Cell *cells, NodeIntType Ncells);

    void saveTree(std::ostream& o) const;
#endif
  protected:
    Node *nodes;
    NodeIntType numNodes;
    NodeIntType lastNode;

    Node *root;
    Cell **sortingHelper;
    Cell *base_cell;

    bool periodic;
    coords replicate;

    Node *buildTree(Cell **cell0,
		    NodeIntType NumCells,
		    uint32_t depth,
		    coords minBound,
		    coords maxBound);
    
    template<bool justCount>
    void recursiveIntersectionCells(RecursionInfoCells<N,ValType, CType>& info,
				    Node *node,
				    int level);

    CoordType computeDistance(const Cell *cell, const coords& x) const;
    void recursiveNearest(Node *node,
			  int level,
			  const coords& x,
			  CoordType& R2,
			  Cell*& cell);
    void recursiveMultipleNearest(RecursionMultipleInfo<N,ValType,CType>& info, Node *node,
				  int level);    

  };

  template<int N, typename T, typename CType>
  NodeIntType gatherActiveCells(KDCell<N,T,CType> **cells, NodeIntType numCells);

};

#include "mykdtree.tcc"

#endif
