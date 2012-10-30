#include <cstring>
#include <algorithm>
#include <limits>
#include <iostream>
#include <cassert>

namespace CosmoTool {

  template<int N, typename ValType, typename CType>
  class CellCompare
  {
  public:
    CellCompare(int k)
    {
      rank = k;
    }

    bool operator()(const KDLeafCell<N,ValType,CType> *a, const KDLeafCell<N,ValType,CType> *b) const
    {
      return (a->coord[rank] < b->coord[rank]);
    }
  protected:
    int rank;
  };

  template<int N, typename ValType, typename CType>
  KDLeafTree<N,ValType,CType>::~KDLeafTree()
  {
  }

  template<int N, typename ValType, typename CType>
  KDLeafTree<N,ValType,CType>::KDLeafTree(Cell *cells, uint32_t Ncells)
  {
    numNodes = Ncells*3;
    numCells = Ncells;
    nodes = new Node[numNodes];

    sortingHelper = new Cell *[Ncells];
    for (uint32_t i = 0; i < Ncells; i++)
	sortingHelper[i] = &cells[i];    

    optimize();
  }

  template<int N, typename ValType, typename CType>
  void KDLeafTree<N,ValType,CType>::optimize()
  {
    coords absoluteMin, absoluteMax;

    std::cout << "Optimizing the tree..." << std::endl;
    uint32_t activeCells = gatherActiveCells(sortingHelper, numCells);
    std::cout << "  number of active cells = " << activeCells << std::endl;

    lastNode = 0;
    for (int i = 0; i < N; i++)
      {
	absoluteMin[i] = std::numeric_limits<typeof (absoluteMin[0])>::max();
	absoluteMax[i] = -std::numeric_limits<typeof (absoluteMax[0])>::max();
      }
    // Find min and max corner
    for (uint32_t i = 0; i < activeCells; i++)
      {
        KDLeafCell<N,ValType,CType> *cell = sortingHelper[i];

        for (int k = 0; k < N; k++) {
          if (cell->coord[k] < absoluteMin[k])
            absoluteMin[k] = cell->coord[k];
          if (cell->coord[k] > absoluteMax[k])
            absoluteMax[k] = cell->coord[k];
        }
      }
    
    std::cout << "   rebuilding the tree..." << std::endl;
    root = buildTree(sortingHelper, activeCells, 0, absoluteMin, absoluteMax);    
    std::cout << "   done." << std::endl;
  }

  template<int N, typename ValType, typename CType>
  uint32_t gatherActiveCells(KDLeafCell<N,ValType,CType> **cells,
			     uint32_t Ncells)
  {
    uint32_t swapId = Ncells-1;
    uint32_t i = 0;

    while (!cells[swapId]->active && swapId > 0)
      swapId--;

    while (i < swapId)
      {
	if (!cells[i]->active)
	  {
	    std::swap(cells[i], cells[swapId]);
	    while (!cells[swapId]->active && swapId > i)
	      {
		swapId--;
	      }
	  }
	i++;
      }
    return swapId+1;
  }

  template<int N, typename ValType, typename CType>
  KDLeafTreeNode<N,ValType,CType> *
        KDLeafTree<N,ValType,CType>::buildTree(Cell **cell0,
					       uint32_t Ncells,
					       uint32_t depth,
					       coords minBound,
					       coords maxBound)
  {
    if (Ncells == 0)
      return 0;

    int axis = depth % N;
    assert(lastNode != numNodes);
    Node *node = &nodes[lastNode++];
    uint32_t mid = Ncells/2;
    coords tmpBound;
    
    // Isolate the environment
    {
      CellCompare<N,ValType,CType> compare(axis);
      std::sort(cell0, cell0+Ncells, compare);
    }

    node->leaf = false;
    memcpy(&node->minBound[0], &minBound[0], sizeof(coords));
    memcpy(&node->maxBound[0], &maxBound[0], sizeof(coords));

    if (Ncells == 1)
      {
	node->leaf = true;
	node->value = *cell0;

#ifdef __KDLEAF_TREE_NUMNODES
	node->numNodes = 1;
#endif
	return node;
      }

    memcpy(tmpBound, maxBound, sizeof(coords));
    tmpBound[axis] = (*(cell0+mid))->coord[axis];
    depth++;
    node->children[0] = buildTree(cell0, mid, depth, minBound, tmpBound);
    
    memcpy(tmpBound, minBound, sizeof(coords));
    tmpBound[axis] = (*(cell0+mid))->coord[axis];
    node->children[1] = buildTree(cell0+mid, Ncells-mid, depth, 
				  tmpBound, maxBound);

#ifdef __KDLEAF_TREE_NUMNODES
    node->numNodes = (node->children[0] != 0) ? node->children[0]->numNodes : 0;
    node->numNodes += (node->children[1] != 0) ? node->children[1]->numNodes : 0;
#endif

    return node;
  }

  template<int N, typename ValType, typename CType>
  uint32_t KDLeafTree<N,ValType,CType>::countActives() const
  {
    uint32_t numActive = 0;
    for (uint32_t i = 0; i < lastNode; i++)
      {
	if (nodes[i].value->active)
	  numActive++;
      }
    return numActive;
  }

  template<int N, typename ValType, typename CType>
  typename KDLeafDef<N,CType>::CoordType
  KDLeafTree<N,ValType,CType>::computeDistance(const Cell *cell, const coords& x) const
  {
    CoordType d2 = 0;

    for (int i = 0; i < N; i++)
      {
	CoordType delta = cell->coord[i] - x[i];
	d2 += delta*delta;
      }
    return d2;
  }

 template<int N, typename ValType, typename CType>
 double KDLeafTree<N,ValType,CType>::countInRange(CType sLo, CType sHigh, Node *root1, Node *root2) const
 {
   double result = recursiveCountInRange((root1 == 0) ? root : root1,
				(root2 == 0) ? root : root2, 
				sLo*sLo, sHigh*sHigh); 
   return result;
 }

 
 template<int N, typename ValType, typename CType>
 double KDLeafTree<N,ValType,CType>::recursiveCountInRange(Node *na, Node *nb,
                                                         CType sLo, CType sHi) const
 {
   assert(nb != 0);
   if (na == 0)
     {
       return 0;
     }
   
   uint32_t numNa = getNumberInNode(na);
   uint32_t numNb = getNumberInNode(nb);
   double Cleft, Cright;
   CType minDist, maxDist;

   if (numNa == 1 && numNb == 1)
     {
       assert(na->leaf && nb->leaf);
       CType ab_dist = computeDistance(na->value, nb->value->coord);
       if (ab_dist >= sLo && ab_dist < sHi)
	 return 1;
       else
	 return 0;
     }
   assert(numNa > 1 || numNb > 1);
   
   bool overlapping_a = true, overlapping_b = true;
   for (int k = 0; k < N; k++) 
     {
       bool min_a_in_B = 
	 ((na->minBound[k] >= nb->minBound[k] && 
	   na->minBound[k] <= nb->maxBound[k]));
       bool max_a_in_B =
	 ((na->maxBound[k] >= nb->minBound[k] && 
	   na->maxBound[k] <= nb->maxBound[k]));
       bool min_b_in_A = 
	 ((nb->minBound[k] >= na->minBound[k] && 
	   nb->minBound[k] <= na->maxBound[k]));
       bool max_b_in_A = 
	 ((nb->maxBound[k] >= na->minBound[k] && 
	   nb->maxBound[k] <= na->maxBound[k]));

      if (!min_a_in_B && !max_a_in_B)
         overlapping_a = false;
      if (!min_b_in_A && !max_b_in_A)
         overlapping_b = false;
     }

   if (overlapping_a || overlapping_b)
     {
        minDist = 0;
        maxDist = 0;
        for (int k = 0; k < N; k++)
        {
          CType delta = max(nb->maxBound[k]-na->minBound[k],na->maxBound[k]-nb->minBound[k]);
          maxDist += delta*delta;
        }
     }
   else
     {
       minDist = maxDist = 0;
       for (int k = 0; k < N; k++)
         {
           CType delta2;
           delta2 = max(nb->maxBound[k]-na->minBound[k],
			na->maxBound[k]-nb->minBound[k]);
           maxDist += delta2*delta2;
         }
       // mins and maxs
       CType minmax[N][2];
       for (int k = 0; k < N; k++)
	 {
	   if (na->minBound[k] < nb->minBound[k])
	     {
	       minmax[k][1] = na->maxBound[k];
	       minmax[k][0] = nb->minBound[k];
	     }
	   else
	     {
	       minmax[k][1] = nb->maxBound[k];
	       minmax[k][0] = na->minBound[k];
	     }
	 }
       for (int k = 0; k < N; k++)
	 {
	   CType delta = max(minmax[k][0]-minmax[k][1], 0.);
	   minDist += delta*delta;
	 }
     }
   
   if (minDist >= sHi)
     return 0;
   if (maxDist < sLo)
     return 0;
   
   if (sLo <= minDist && maxDist < sHi)
     return ((double)numNa)*numNb;
   
   if (numNa < numNb) 
     {
       assert(!nb->leaf);
       Cleft = recursiveCountInRange(nb->children[0], na, sLo, sHi);
       Cright = recursiveCountInRange(nb->children[1], na, sLo, sHi);
     }
   else
     {
       assert(!na->leaf);
       Cleft = recursiveCountInRange(na->children[0], nb, sLo, sHi);
       Cright = recursiveCountInRange(na->children[1], nb, sLo, sHi);
     }
   return Cleft+Cright;
 }

};
