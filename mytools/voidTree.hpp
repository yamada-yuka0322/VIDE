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

  int lookupParent(int voidId)
  {
    std::set<ZobovZone *> sref;
    sref.insert(zobov.allVoids[voidId].zId.begin(), zobov.allVoids[voidId].zId.end());
    
    
    std::vector<ZobovZone *> sout(sref.size());  
    int lastSize = 0x7fffffff;
    int goodParent = voidId;

    // Voids are sorted according to volume
    for (int i = voidId-1; i >= 0; i--) 
      {
	std::set<ZobovZone *> s1;
	
	int counter = 0;
	s1.insert(zobov.allVoids[i].zId.begin(), zobov.allVoids[i].zId.end());
	for (std::set<ZobovZone *>::iterator iter = s1.begin(), iter2 = sref.begin();
	     iter != s1.end() && iter2 != sref.end();
	     ++iter)
	  {
	    if (*iter == *iter2)
	      {
		counter++;
		++iter2;
	      }
	    else while (*iter > *iter2 && iter2 != sref.end())
		   ++iter2;
	  }
	
	if (counter == sref.size() && s1.size() < lastSize)
	  {
	    return i;
	  }
      }
    return -1;
  }
  
  VoidTree(ZobovRep& rep)
    : zobov(rep)
  {
    totalNumNodes = rep.allVoids.size();
    
    nodes = new VoidNode[totalNumNodes];

    for (int i = 0; i < rep.allVoids.size(); i++)
      {
	nodes[i].vid = i;
	nodes[i].parent = 0;
	nodes[i].children.clear();
      }

    std::cout << "Linking voids together..." << std::endl;
    double volMin = 4*M_PI/3*pow(4.*512/500.,3);
    int inserted = 0;
    for (int i = rep.allVoids.size()-1; i>=0;i--)
      {
        if (rep.allVoids[i].volume < volMin) continue;

	int p = lookupParent(i);

	std::cout << i << std::endl;
	
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
    for (int i = inserted; i < totalNumNodes; i++)
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
    VoidList::iterator i;
    
    if (!traverse(node))
      return;

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
