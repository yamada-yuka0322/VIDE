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
