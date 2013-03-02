/*+
    VIDE -- Void IDEntification pipeline -- ./c_tools/libzobov/voidTree.hpp
    Copyright (C) 2010-2013 Guilhem Lavaux
    Copyright (C) 2011-2013 Paul M. Sutter

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; version 2 of the License.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
+*/


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

struct VoidNodeOnDisk
{
  int vid;
  int parent;
};

class VoidTree
{
protected:
  uint32_t totalNumNodes, activeNodes;
  VoidNode *nodes;
  VoidNode *rootNode;  
  ZobovRep& zobov;
public:
  typedef std::list<VoidNode *> VoidList;

  void dumpTree(std::ostream& o)
  {
    VoidNodeOnDisk data;

    o.write((char*)&activeNodes, sizeof(uint32_t));
    for (uint32_t i = 0; i < activeNodes; i++)
      {
	data.vid = nodes[i].vid;
	if (nodes[i].parent == 0)
	  data.parent = -1;
	else
	  data.parent = nodes[i].parent - nodes;
	o.write((char *)&data, sizeof(VoidNodeOnDisk));
      }
    
  }  

  int lookupParent(int voidId, const std::vector<std::list<int> >& voids_for_zones)
  {
    int lastSize = 0x7fffffff;
    int goodParent = -1;
    ZobovVoid &ref_void = zobov.allVoids[voidId];
    const std::list<int>& candidateList = voids_for_zones[ref_void.zId.front()];
    std::list<int>::const_iterator iter_candidate = candidateList.begin();

    //    std::cout << "candidate list size is " << candidateList.size() << std::endl;

    while (iter_candidate != candidateList.end())
      {
        int vid_candidate = *iter_candidate;

        if (vid_candidate == voidId)
          break;
        ++iter_candidate;
      }
    if (iter_candidate == candidateList.end())
      {
//        std::cout << "Failure to lookup parent" << std::endl;
        return -1;
      }

    // voidId must be in the list.
    //    assert(iter_candidate != candidateList.end());
    // Go back

    iter_candidate = candidateList.end();

    int vid_good_candidate = -1;
    int old_good_candidate_size = zobov.allZones.size()+1;

    do
      {
        int vid_candidate;

        --iter_candidate;

        vid_candidate = *iter_candidate;
        std::vector<int>& candidate_zIds = zobov.allVoids[vid_candidate].zId;

	if (voidId == vid_candidate)
	  continue;

        if (candidate_zIds.size() < ref_void.zId.size())
          {
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
	  {
	    if (candidate_zIds.size() < old_good_candidate_size)
	      {
		vid_good_candidate = vid_candidate;
		old_good_candidate_size = candidate_zIds.size();
	      }
	    //	    std::cout << "Found parent " << vid_candidate << std::endl;
	    //	    return vid_candidate;
	  }
	
        // Go bigger, though I would say we should not to.
      }
    while (iter_candidate != candidateList.begin()) ;
    //if (vid_good_candidate < 0)
    //  std::cout << "Failure to lookup parent (2)" << std::endl;
    return vid_good_candidate;
  }

  VoidTree(ZobovRep& rep, std::istream& disk)
    : zobov(rep)
  {
    totalNumNodes = rep.allVoids.size();

    disk.read((char *)&activeNodes, sizeof(uint32_t));
    nodes = new VoidNode[activeNodes];
    rootNode = 0;
    for (uint32_t i = 0; i < activeNodes; i++)
      {
	VoidNodeOnDisk data;

	disk.read((char *)&data, sizeof(data));	
	nodes[i].vid = data.vid;
	if (data.parent < 0)
	  {
	    if (rootNode != 0)
	      {
		std::cerr << "Multiple root to the tree !!!" << std::endl;
		abort();
	      }
	    nodes[i].parent = 0;
	    rootNode = &nodes[i];
	  }
	else
	  {
	    nodes[i].parent = nodes + data.parent;
	    nodes[i].parent->children.push_back(&nodes[i]);
	  }
      }
    
    computeMaxDepth();
    computeChildrenByNode();
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

    // One additional for the mega-root
    nodes = new VoidNode[totalNumNodes+1];

    for (int i = 0; i <= totalNumNodes; i++)
      {
	nodes[i].vid = i;
	nodes[i].parent = 0;
      }

    std::cout << "Linking voids together..." << std::endl;
    double volMin = 0;// 4*M_PI/3*pow(4.*512/500.,3);
    int inserted = 0;
    for (int i = 0; i < rep.allVoids.size(); i++)
      {
        if (rep.allVoids[i].volume < volMin) continue;

	int p = lookupParent(i, voids_for_zones);
        if ((i % 1000) == 0) std::cout << i << std::endl;

	if (p >= 0)
	  {
	    nodes[p].children.push_back(&nodes[i]);
	    nodes[i].parent = &nodes[p];
	  }
	  
	inserted++;
      }

    assert(inserted <= totalNumNodes); 
    rootNode = &nodes[inserted];
    rootNode->vid = -1;
    rootNode->parent = 0;

    for (int i = 0; i < inserted; i++)
      if (nodes[i].parent == 0)
	{
	  nodes[i].parent = rootNode;
	  rootNode->children.push_back(&nodes[i]);
	}
    activeNodes = inserted+1;
    
    computeMaxDepth();
    computeChildrenByNode();
  }

  ~VoidTree()
  {
    delete[] nodes;
  }

  int _depth_computer(VoidNode *node)
  {
    VoidList::iterator i = node->children.begin();
    int d = 0;

    while (i != node->children.end())
      {
	d = std::max(d,_depth_computer(*i));
	++i;
      }  

    return d+1;
  }

  void computeMaxDepth()
  {
    std::cout << "maximum depth is " <<  _depth_computer(rootNode) << std::endl;
  }

  struct _children_stat {
    int num, min_num, max_num, num_zero,num_one, num_multiple;
  };

  void _children_computer(VoidNode *node, _children_stat& s)
  {
    VoidList::iterator i = node->children.begin();
    int d = 0, j = 0;

    while (i != node->children.end())
      {
        _children_computer(*i, s);
        ++i;
        ++j;
      }
    s.num += j;
    if (j!= 0)
      s.min_num = std::min(s.min_num, j);
    else s.num_zero ++;
    if (j==1) s.num_one++;
    if (j>1) s.num_multiple++;
    s.max_num = std::max(s.max_num, j);
  }

  void computeChildrenByNode()
  {
    _children_stat s;
   s.num = 0;
   s.min_num = activeNodes+1;
   s.max_num = s.num_zero = s.num_one =s.num_multiple= 0;
   _children_computer(rootNode, s);
    std::cout << "Average children by node " << s.num*1.0/activeNodes << " , " << s.min_num << " " << s.max_num << " " << s.num_zero << " " << s.num_one << " " << s.num_multiple << std::endl;
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


  VoidNode *getRoot() { return rootNode; }

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

  template<typename T, typename T2>
  void walkNodeWithMark(VoidNode *node, T& traverse, const T2& mark)
  {
    T2 new_mark = mark;

    if (!traverse(node, new_mark))
      return;

    VoidList::iterator i = node->children.begin();
    
    while (i != node->children.end())
      {
	walkNodeWithMark(*i, traverse, new_mark);
	++i;
      }    
  }

  template<typename T,typename T2>
  void walkWithMark(T& traverse, T2 mark)
  {
    walkNodeWithMark(rootNode, traverse, mark);
  }
};

#endif
