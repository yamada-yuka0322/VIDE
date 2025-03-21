#include "replicateGenerator.hpp"
#include <cstring>
#include "omptl/omptl_algorithm"
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

    bool operator()(const KDCell<N,ValType,CType> *a, const KDCell<N,ValType,CType> *b) const
    {
      return (a->coord[rank] < b->coord[rank]);
    }
  protected:
    int rank;
  };

  template<int N, typename ValType, typename CType, typename CellSplitter>
  KDTree<N,ValType,CType,CellSplitter>::~KDTree()
  {
  }

  template<int N, typename ValType, typename CType, typename CellSplitter>
  KDTree<N,ValType,CType,CellSplitter>::KDTree(Cell *cells, NodeIntType Ncells)
  {
    periodic = false;

    base_cell = cells;
    numNodes = Ncells;
    nodes = new Node[numNodes];

    sortingHelper = new Cell *[Ncells];
    for (NodeIntType i = 0; i < Ncells; i++)
	sortingHelper[i] = &cells[i];

    optimize();
  }

  template<int N, typename ValType, typename CType, typename CellSplitter>
  void KDTree<N,ValType,CType,CellSplitter>::optimize()
  {
    coords absoluteMin, absoluteMax;

    std::cout << "Optimizing the tree..." << std::endl;
    NodeIntType activeCells = gatherActiveCells(sortingHelper, numNodes);
    std::cout << "  number of active cells = " << activeCells << std::endl;

    lastNode = 0;
    for (int i = 0; i < N; i++)
      {
	absoluteMin[i] = std::numeric_limits<typeof (absoluteMin[0])>::max();
	absoluteMax[i] = -std::numeric_limits<typeof (absoluteMax[0])>::max();
      }
    // Find min and max corner
    for (NodeIntType i = 0; i < activeCells; i++)
      {
        KDCell<N,ValType,CType> *cell = sortingHelper[i];

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

  template<int N, typename ValType, typename CType, typename CellSplitter>
  uint64_t KDTree<N,ValType,CType,CellSplitter>::getIntersection(const coords& x, CoordType r,
						    KDTree<N,ValType,CType,CellSplitter>::Cell **cells,
						    uint64_t numCells)
  {
    RecursionInfoCells<N,ValType,CType> info;

    memcpy(info.x, x, sizeof(x));
    info.r = r;
    info.r2 = r*r;
    info.cells = cells;
    info.currentRank = 0;
    info.numCells = numCells;
    info.distances = 0;

    recursiveIntersectionCells<false>(info, root, 0);
    if (periodic)
      {
        ReplicateGenerator<float, N> r(x, replicate);

        do
          {
            coords x_new;
            r.getPosition(info.x);
            recursiveIntersectionCells<false>(info, root, 0);
          }
        while (r.next());
      }

    return info.currentRank;
  }

  template<int N, typename ValType, typename CType, typename CellSplitter>
  uint64_t KDTree<N,ValType,CType,CellSplitter>::getIntersection(const coords& x, CoordType r,
						    Cell **cells,
						    CoordType *distances,
						    uint64_t numCells)
  {
    RecursionInfoCells<N,ValType> info;

    memcpy(info.x, x, sizeof(x));
    info.r = r;
    info.r2 = r*r;
    info.cells = cells;
    info.currentRank = 0;
    info.numCells = numCells;
    info.distances = distances;

    recursiveIntersectionCells<false>(info, root, 0);
    if (periodic)
      {
        ReplicateGenerator<float, N> r(x, replicate);

        do
          {
            coords x_new;
            r.getPosition(info.x);
            recursiveIntersectionCells<false>(info, root, 0);
          }
        while (r.next());
      }
    return info.currentRank;
  }

  template<int N, typename ValType, typename CType, typename CellSplitter>
  uint64_t KDTree<N,ValType,CType,CellSplitter>::countCells(const coords& x, CoordType r)
  {
    RecursionInfoCells<N,ValType> info;

    memcpy(info.x, x, sizeof(x));
    info.r = r;
    info.r2 = r*r;
    info.cells = 0;
    info.currentRank = 0;
    info.numCells = 0;
    info.distances = 0;

    recursiveIntersectionCells<true>(info, root, 0);
    if (periodic)
      {
        ReplicateGenerator<float, N> r(x, replicate);

        do
          {
            coords x_new;
            r.getPosition(info.x);
            recursiveIntersectionCells<true>(info, root, 0);
          }
        while (r.next());
      }

    return info.currentRank;
  }

  template<int N, typename ValType, typename CType, typename CellSplitter>
  template<bool justCount>
  void KDTree<N,ValType,CType,CellSplitter>::recursiveIntersectionCells(RecursionInfoCells<N,ValType,CType>& info,
							   Node *node,
							   int level)
  {
    int axis = level % N;
    CoordType d2 = 0;

#if __KD_TREE_ACTIVE_CELLS == 1
    if (node->value->active)
#endif
      {
	for (int j = 0; j < 3; j++)
	  {
	    CoordType delta = info.x[j]-node->value->coord[j];
	    d2 += delta*delta;
	  }
	if (d2 < info.r2)
	  {
	    if (!justCount)
	      {
		if (info.currentRank == info.numCells)
		  throw NotEnoughCells();
		info.cells[info.currentRank] = node->value;
		if (info.distances)
		  info.distances[info.currentRank] = d2;
	      }
	    info.currentRank++;
	  }
      }

    // The hypersphere intersects the left child node
    if (((info.x[axis]+info.r) > node->minBound[axis]) &&
	((info.x[axis]-info.r) < node->value->coord[axis]))
      {
	if (node->children[0] != 0)
	  recursiveIntersectionCells<justCount>(info, node->children[0],
				     level+1);
      }
    if (((info.x[axis]+info.r) > node->value->coord[axis]) &&
	((info.x[axis]-info.r) < node->maxBound[axis]))
      {
	if (node->children[1] != 0)
	  recursiveIntersectionCells<justCount>(info, node->children[1],
					       level+1);
      }
  }

  template<int N, typename ValType, typename CType>
  NodeIntType gatherActiveCells(KDCell<N,ValType,CType> **cells,
			        NodeIntType Ncells)
  {
    NodeIntType swapId = Ncells-1;
    NodeIntType i = 0;

#if __KD_TREE_ACTIVE_CELLS == 1
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
#endif
    return swapId+1;
  }

  template<int N, typename ValType, typename CType>
  void KD_default_cell_splitter<N,ValType,CType>::operator()(KDCell<N,ValType,CType> **cells, NodeIntType Ncells,
                                                             NodeIntType& split_index, int axis,
                                                             typename KDDef<N,CType>::KDCoordinates minBound,
                                                             typename KDDef<N,CType>::KDCoordinates maxBound)
  {
    CellCompare<N,ValType,CType> compare(axis);
    omptl::sort(cells,cells+Ncells,compare); // std::sort(cells, cells+Ncells, compare);
    split_index = Ncells/2;
  }


  template<int N, typename ValType, typename CType, typename CellSplitter>
  KDTreeNode<N,ValType,CType> *KDTree<N,ValType,CType,CellSplitter>::buildTree(Cell **cell0,
								  NodeIntType Ncells,
								  uint32_t depth,
								  coords minBound,
								  coords maxBound)
  {
    if (Ncells == 0)
      return 0;

    Node *node;
    int axis = depth % N;
    NodeIntType mid;
    coords tmpBound;
    NodeIntType nodeId;

//#pragma omp atomic capture
    nodeId = (this->lastNode)++;

    node = &nodes[nodeId];

    // Isolate the environment
    splitter(cell0, Ncells, mid, axis, minBound, maxBound);

    node->value = *(cell0+mid);
    memcpy(&node->minBound[0], &minBound[0], sizeof(coords));
    memcpy(&node->maxBound[0], &maxBound[0], sizeof(coords));

    memcpy(tmpBound, maxBound, sizeof(coords));
    tmpBound[axis] = node->value->coord[axis];

    depth++;
//#pragma omp task private(tmpBound)
    {
      node->children[0] = buildTree(cell0, mid, depth, minBound, tmpBound);
    }

    memcpy(tmpBound, minBound, sizeof(coords));
    tmpBound[axis] = node->value->coord[axis];
#pragma omp task private(tmpBound)
    {
      node->children[1] = buildTree(cell0+mid+1, Ncells-mid-1, depth,
				  tmpBound, maxBound);
    }

#pragma omp taskwait

#ifdef __KD_TREE_NUMNODES
    node->numNodes = (node->children[0] != 0) ? node->children[0]->numNodes : 0;
    node->numNodes += (node->children[1] != 0) ? node->children[1]->numNodes : 0;
    node->numNodes++;
#endif

    return node;
  }

  template<int N, typename ValType, typename CType, typename CellSplitter>
  NodeIntType KDTree<N,ValType,CType,CellSplitter>::countActives() const
  {
    NodeIntType numActive = 0;
    for (NodeIntType i = 0; i < lastNode; i++)
      {
#if __KD_TREE_ACTIVE_CELLS == 1
	if (nodes[i].value->active)
#endif
	  numActive++;
      }
    return numActive;
  }

  template<int N, typename ValType, typename CType, typename CellSplitter>
  typename KDDef<N,CType>::CoordType
  KDTree<N,ValType,CType,CellSplitter>::computeDistance(const Cell *cell, const coords& x) const
  {
    CoordType d2 = 0;

    for (int i = 0; i < N; i++)
      {
	CoordType delta = cell->coord[i] - x[i];
	d2 += delta*delta;
      }
    return d2;
  }

  template<int N, typename ValType, typename CType, typename CellSplitter>
  void
  KDTree<N,ValType,CType,CellSplitter>::recursiveNearest(
				       Node *node,
				       int level,
				       const coords& x,
				       CoordType& R2,
				       Cell *& best)
  {
    CoordType d2 = 0;
    int axis = level % N;
    Node *other, *go;

    if (x[axis] < node->value->coord[axis])
      {
	// The best is potentially in 0.
	go = node->children[0];
	other = node->children[1];
      }
    else
      {
	// If not it is in 1.
	go = node->children[1];
	other = node->children[0];
	if (go == 0)
	  {
	    go = other;
	    other = 0;
	  }
      }

    if (go != 0)
      {
	recursiveNearest(go, level+1,
			 x, R2,best);
      }
    else
      {
	CoordType thisR2 = computeDistance(node->value, x);
	if (thisR2 < R2)
	  {
	    R2 = thisR2;
	    best = node->value;
	  }
	return;
      }

    // Check if current node is not the nearest
    CoordType thisR2 =
      computeDistance(node->value, x);

    if (thisR2 < R2)
      {
	R2 = thisR2;
	best = node->value;
      }

    // Now we found the best. We check whether the hypersphere
    // intersect the hyperplane of the other branch

    CoordType delta1;

    delta1 = x[axis]-node->value->coord[axis];
    if (delta1*delta1 < R2)
      {
	// The hypersphere intersects the hyperplane. Try the
	// other branch
	if (other != 0)
	  {
	    recursiveNearest(other, level+1, x, R2, best);
	  }
      }
  }

 template<int N, typename ValType, typename CType, typename CellSplitter>
 KDCell<N,ValType,CType> *
 KDTree<N,ValType,CType,CellSplitter>::getNearestNeighbour(const coords& x)
 {
   CoordType R2 = INFINITY;
   Cell *best = 0;

   recursiveNearest(root, 0, x, R2, best);
   if (periodic)
     {
       ReplicateGenerator<float, N> r(x, replicate);

       do
         {
           coords x_new;
           r.getPosition(x_new);
           recursiveNearest(root, 0, x_new, R2, best);
         }
       while (r.next());
     }

   return best;
 }

  template<int N, typename ValType, typename CType, typename CellSplitter>
  void
  KDTree<N,ValType,CType,CellSplitter>::recursiveMultipleNearest(RecursionMultipleInfo<N,ValType,CType>& info, Node *node,
					       int level)
  {
    CoordType d2 = 0;
    int axis = level % N;
    Node *other, *go;

    if (info.x[axis] < node->value->coord[axis])
      {
	// The best is potentially in 0.
	go = node->children[0];
	other = node->children[1];
      }
    else
      {
	// If not it is in 1.
	go = node->children[1];
	other = node->children[0];
	//	if (go == 0)
	//	  {
	//	go = other;
	//other = 0;
	//}
      }

    if (go != 0)
      {
	recursiveMultipleNearest(info, go, level+1);
      }

    // Check if current node is not the nearest
    CoordType thisR2 =
      computeDistance(node->value, info.x);
    info.queue.push(node->value, thisR2);
    info.traversed++;
    //    if (go == 0)
    //      return;

    // Now we found the best. We check whether the hypersphere
    // intersect the hyperplane of the other branch

    CoordType delta1;

    delta1 = info.x[axis]-node->value->coord[axis];
    if (delta1*delta1 < info.queue.getMaxPriority())
      {
	// The hypersphere intersects the hyperplane. Try the
	// other branch
	if (other != 0)
	  {
	    recursiveMultipleNearest(info, other, level+1);
	  }
      }
  }

 template<int N, typename ValType, typename CType, typename CellSplitter>
 void KDTree<N,ValType,CType,CellSplitter>::getNearestNeighbours(const coords& x, uint64_t N2,
						    Cell **cells)
 {
   RecursionMultipleInfo<N,ValType> info(x, cells, N2);

   for (int i = 0; i < N2; i++)
     cells[i] = 0;

   recursiveMultipleNearest(info, root, 0);
   if (periodic)
    {
       ReplicateGenerator<float, N> r(x, replicate);

       do
         {
           coords x_new;
           r.getPosition(info.x);
           recursiveMultipleNearest(info, root, 0);
         }
       while (r.next());
    }

   //   std::cout << "Traversed = " << info.traversed << std::endl;
 }

 template<int N, typename ValType, typename CType, typename CellSplitter>
 void KDTree<N,ValType,CType,CellSplitter>::getNearestNeighbours(const coords& x, uint64_t N2,
						    Cell **cells,
						    CoordType *distances)
 {
   RecursionMultipleInfo<N,ValType> info(x, cells, N2);

   for (int i = 0; i < N2; i++)
     cells[i] = 0;

   recursiveMultipleNearest(info, root, 0);
   if (periodic)
    {
       ReplicateGenerator<float, N> r(x, replicate);

       do
         {
           coords x_new;
           r.getPosition(info.x);
           recursiveMultipleNearest(info, root, 0);
         }
       while (r.next());
    }
   memcpy(distances, info.queue.getPriorities(), sizeof(CoordType)*N2);
 }

#ifdef __KD_TREE_SAVE_ON_DISK
#define KDTREE_DISK_SIGNATURE "KDTREE"
#define KDTREE_DISK_SIGNATURE_LEN 7

  template<int N, typename CType>
  struct KDTreeOnDisk
  {
    long cell_id;
    long children_node[2];
    typename KDDef<N, CType>::KDCoordinates minBound, maxBound;
  };

  struct KDTreeHeader
  {
    char id[KDTREE_DISK_SIGNATURE_LEN];
    long nodesUsed, numCells;
    long rootId;
  };

  template<int N, typename ValType, typename CType, typename CellSplitter>
  void KDTree<N,ValType,CType,CellSplitter>::saveTree(std::ostream& o) const
  {
    KDTreeHeader h;

    strncpy(h.id, KDTREE_DISK_SIGNATURE, KDTREE_DISK_SIGNATURE_LEN);
    h.nodesUsed = lastNode;
    h.numCells = numNodes;
    h.rootId = root - nodes;
    o.write((char*)&h, sizeof(h));

    for (long i = 0; i < lastNode; i++)
      {
	KDTreeOnDisk<N,CType> node_on_disk;

	node_on_disk.cell_id = nodes[i].value - base_cell;
	if (nodes[i].children[0] == 0)
	  node_on_disk.children_node[0] = -1L;
	else
	  node_on_disk.children_node[0] = nodes[i].children[0] - nodes;
       assert((node_on_disk.children_node[0] == -1) || ((node_on_disk.children_node[0] >= 0) &&  (node_on_disk.children_node[0] < lastNode)));

	if (nodes[i].children[1] == 0)
	  node_on_disk.children_node[1] = -1L;
	else
	  node_on_disk.children_node[1] = nodes[i].children[1] - nodes;
       assert((node_on_disk.children_node[1] == -1) || ((node_on_disk.children_node[1] >= 0) &&  (node_on_disk.children_node[1] < lastNode)));

	memcpy(node_on_disk.minBound, nodes[i].minBound, sizeof(coords));
	memcpy(node_on_disk.maxBound, nodes[i].maxBound, sizeof(coords));

	o.write((char *)&node_on_disk, sizeof(node_on_disk));
      }
  }

  template<int N, typename ValType, typename CType, typename CellSplitter>
  KDTree<N,ValType,CType,CellSplitter>::KDTree(std::istream& in, Cell *cells, NodeIntType Ncells)
  {
    KDTreeHeader h;

    if (!in)
      throw InvalidOnDiskKDTree();

    in.read((char *)&h, sizeof(h));
    if (!in || strncmp(h.id, KDTREE_DISK_SIGNATURE, KDTREE_DISK_SIGNATURE_LEN) != 0)
      {
         std::cerr << "KDTree Signature invalid" << std::endl;
         throw InvalidOnDiskKDTree();
      }

    if (h.numCells != Ncells || h.nodesUsed < 0) {
      std::cerr << "The number of cells has changed (" << h.numCells << " != " << Ncells << ") or nodesUsed=" << h.nodesUsed << std::endl;
      throw InvalidOnDiskKDTree();
    }

    base_cell = cells;
    nodes = new Node[h.nodesUsed];
    lastNode = h.nodesUsed;
    numNodes = Ncells;

    for (long i = 0; i < lastNode; i++)
      {
	KDTreeOnDisk<N,CType> node_on_disk;

	in.read((char *)&node_on_disk, sizeof(node_on_disk));

	if (!in) {
           std::cerr << "End-of-file reached" << std::endl;
           delete[] nodes;
           throw InvalidOnDiskKDTree();
        }
        if (node_on_disk.cell_id > numNodes || node_on_disk.cell_id < 0 ||
	    (node_on_disk.children_node[0] >= 0 && node_on_disk.children_node[0] > lastNode) || node_on_disk.children_node[0] < -1 ||
	    (node_on_disk.children_node[1] >= 0 && node_on_disk.children_node[1] > lastNode) || node_on_disk.children_node[1] < -1)
	  {
	    delete[] nodes;
            std::cerr << "Invalid cell id or children node id invalid" << std::endl;
            std::cerr << node_on_disk.cell_id << std::endl << node_on_disk.children_node[0] << std::endl << node_on_disk.children_node[1] << std::endl;
	    throw InvalidOnDiskKDTree();
	  }

	nodes[i].value = base_cell + node_on_disk.cell_id;
	if (node_on_disk.children_node[0] == -1)
	  nodes[i].children[0] = 0;
	else
	  nodes[i].children[0] = nodes + node_on_disk.children_node[0];

	if (node_on_disk.children_node[1] == -1)
	  nodes[i].children[1] = 0;
	else
	  nodes[i].children[1] = nodes + node_on_disk.children_node[1];

	memcpy(nodes[i].minBound, node_on_disk.minBound, sizeof(coords));
	memcpy(nodes[i].maxBound, node_on_disk.maxBound, sizeof(coords));

	int c;
	for (c = 0; c < N; c++)
	  if (nodes[i].value->coord[c] < nodes[i].minBound[c] ||
	      nodes[i].value->coord[c] > nodes[i].maxBound[c])
	    break;
	if (c != N)
	  {
	    delete[] nodes;
            std::cerr << "Coordinates of the cell inconsistent with the boundaries" << std::endl
                      << "   X=" << nodes[i].value->coord[0] << " B=[" << nodes[i].minBound[0] << "," << nodes[i].maxBound[0] << "]" << std::endl
                      << "   Y=" << nodes[i].value->coord[1] << " B=[" << nodes[i].minBound[1] << "," << nodes[i].maxBound[1] << "]" << std::endl
                      << "   Z=" << nodes[i].value->coord[2] << " B=[" << nodes[i].minBound[2] << "," << nodes[i].maxBound[2] << "]" << std::endl;
	    throw InvalidOnDiskKDTree();
	  }
      }

    root = &nodes[h.rootId];

    sortingHelper = new Cell *[Ncells];
    for (NodeIntType i = 0; i < Ncells; i++)
      sortingHelper[i] = &cells[i];
  }
#endif

};
