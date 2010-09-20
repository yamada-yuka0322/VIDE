#ifndef _VOID_TREE_HPP
#define _VOID_TREE_HPP

#include <iostream>
#include <stdint.h>
#include "loadZobov.hpp"
#include <list>
#include <set>
#include <vector>

struct VoidNode
{
  int vid;
  VoidNode *parent;
  std::list<VoidNode *> children;
};

class VoidTree
{
protected:
  uint32_t totalNumNodes;
  VoidNode *nodes;
  VoidNode *rootNode;  
  ZobovRep& zobov;
public:
  typedef std::list<VoidNode *> VoidList;

  int lookupParent(int voidId, const std::vector<std::list<int> >& voids_for_zones)
  {
    int lastSize = 0x7fffffff;
    int goodParent = -1;
    ZobovVoid &ref_void = zobov.allVoids[voidId];
    const std::list<int>& candidateList = voids_for_zones[ref_void.zId.front()];
    std::list<int>::const_iterator iter_candidate = candidateList.begin();

    //std::cout << "candidate list size is " << candidateList.size() << std::endl;

    while (iter_candidate != candidateList.end())
      {
        int vid_candidate = *iter_candidate;

        if (vid_candidate == voidId)
          break;
        ++iter_candidate;
      }
    if (iter_candidate == candidateList.begin())
      {
        std::cout << "Failure to lookup parent" << std::endl;
        return -1;
      }

    // voidId must be in the list.
    assert(iter_candidate != candidateList.end());
    // Go back

    do
      {
        int vid_candidate;

        --iter_candidate;

        vid_candidate = *iter_candidate;
        std::vector<int>& candidate_zIds = zobov.allVoids[vid_candidate].zId;

        if (candidate_zIds.size() < ref_void.zId.size())
          {
            --iter_candidate;
            continue;
          }

	int counter = 0;
        // All zones id are sorted in each void. So we just have parse the 
        // vector of zones and check whether all the zones in ref_void.zId is
        // in iter_candidate->zId, the list is analyzed only once.
        // THOUGHT: candidateList may contain directly the information. It would suffice to have the void ids sorted according to volume. Then we just have to jump to the indice just smaller than voidId.

        int k = 0;
        for (int j = 0; j < candidate_zIds.size() && k < ref_void.zId.size(); j++)
          {
            if (candidate_zIds[j] == ref_void.zId[k])
              k++;
            else if (candidate_zIds[j] > ref_void.zId[k])
              break;
          }
        if (k==ref_void.zId.size())
          return vid_candidate;
	
        // Go bigger, though I would say we should not to.
      }
    while (iter_candidate != candidateList.begin()) ;
    std::cout << "Failure to lookup parent (2)" << std::endl;
    return -1;
  }
  
  VoidTree(ZobovRep& rep)
    : zobov(rep)
  {
    totalNumNodes = rep.allVoids.size();

    std::vector<std::list<int> > voids_for_zones;

    voids_for_zones.resize(rep.allZones.size());
    for (int i = 0; i < rep.allVoids.size(); i++) 
       {
         ZobovVoid& v = rep.allVoids[i];
         for (int j = 0; j < v.zId.size(); j++)
           voids_for_zones[v.zId[j]].push_back(i);
       }

    nodes = new VoidNode[totalNumNodes];

    for (int i = 0; i < rep.allVoids.size(); i++)
      {
	nodes[i].vid = i;
	nodes[i].parent = 0;
	nodes[i].children.clear();
      }

    std::cout << "Linking voids together..." << std::endl;
   double volMin = 0;// 4*M_PI/3*pow(4.*512/500.,3);
    int inserted = 0;
    for (int i = rep.allVoids.size()-1; i>=0;i--)
      {
        if (rep.allVoids[i].volume < volMin) continue;

	int p = lookupParent(i, voids_for_zones);
        if ((i % 1000) == 0) std::cout << i << std::endl;

	if (p < 0)
	  {
	    if (i != 0)
	      std::cerr << "Warning ! Voids without parent and it is not the root !" << std::endl;
	    continue;
	  }

	nodes[p].children.push_back(&nodes[i]);
	nodes[i].parent = &nodes[p];
	inserted++;
      }

    rootNode = 0;
    for (int i = 0; i < inserted; i++)
      if (nodes[i].parent == 0)
	{
	  if (rootNode != 0)
	    {
	      std::cerr << "Multiple root to the tree !!!" << std::endl;
	      abort();
	    }
	  rootNode = &nodes[i];
	}
  }

  ~VoidTree()
  {
    delete[] nodes;
  }

  int getParent(int vid) const
  {
    assert(nodes[vid].parent != 0);
    return nodes[vid].parent->vid;
  }
  
  const VoidList& getChildren(int vid) const
  {
    return nodes[vid].children;
  }


  template<typename T>
  void walkNode(VoidNode *node, T& traverse)
  {
    if (!traverse(node))
      return;

    VoidList::iterator i = node->children.begin();
    
    while (i != node->children.end())
      {
	walkNode(*i, traverse);
	++i;
      }
  }

  template<typename T>
  void walk(T& traverse)
  {
    walkNode(rootNode, traverse);
  }
  
};

#endif
