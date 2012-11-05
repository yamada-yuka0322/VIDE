#define __KD_TREE_SAVE_ON_DISK
#include <ctime>
#include <cassert>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include "mykdtree.hpp"

#define NTRY 50000
#define ND 2

using namespace std;
using namespace CosmoTool;

typedef KDTree<ND,char> MyTree;
typedef MyTree::Cell MyCell;

MyCell *findNearest(MyTree::coords& xc, MyCell *cells, uint32_t Ncells)
{
  MyCell *near2 = 0;
  double R2 = INFINITY;
  for (int i = 0; i < Ncells; i++)
    {
      double d2 = 0;
      for (int j = 0; j < ND; j++)
	{
	  double delta = xc[j]-cells[i].coord[j];
	  d2 += delta*delta;
	}
      if (d2 < R2)
	{
	  near2 = &cells[i];
	  near2->val = i;
	  R2 = d2;
	}
    }
  return near2;
}

int main()
{
  uint32_t Ncells = 3000;
  MyCell *cells = new MyCell[Ncells];
  
  for (int i = 0; i < Ncells; i++)
    {
      cells[i].active = true;
      for (int l = 0; l < ND; l++)
	cells[i].coord[l] = drand48();
    }
  
  MyTree tree(cells, Ncells);
  
  MyTree::coords *xc = new MyTree::coords[NTRY];

  cout << "Generating seeds..." << endl;
  for (int k = 0; k < NTRY; k++)
    {
      for (int l = 0; l < ND; l++)
	xc[k][l] = drand48();
    }

  // Check consistency
  cout << "Check consistency..." << endl;
#if 0
  for (int k = 0; k < NTRY; k++) {
    MyCell *near = tree.getNearestNeighbour(xc[k]);
    MyCell *near2 = findNearest(xc[k], cells, Ncells);
    assert(near == near2);
  }
#endif  
  cout << "Check timing..." << endl;
  // Check timing
  clock_t startTimer = clock();
  
  for (int k = 0; k < NTRY; k++) {
    MyCell *near = tree.getNearestNeighbour(xc[k]);
    near->val = 0;
  }

  clock_t endTimer = clock();
  clock_t delta = endTimer-startTimer;
  double myTime = delta*1.0/CLOCKS_PER_SEC * 1000000.0;

  cout << "KDTree search/sec = " << myTime/NTRY << " us" << endl;
  
  startTimer = clock();
  for (int k = 0; k < NTRY; k++) {
    MyCell *near = findNearest(xc[k], cells, Ncells);
  }
  endTimer = clock();
  delta = endTimer-startTimer;
  myTime = delta*1.0/CLOCKS_PER_SEC * 1000000.0;

  cout << "Direct search/sec = " << myTime/NTRY << " us" << endl;

  {
    ofstream f("kd.dat");
    tree.saveTree(f);
  }

  cout << "Trying to reload the tree" << endl;
  {
    ifstream f("kd.dat");
    MyTree tree2(f, cells, Ncells);
    cout << "Check timing..." << endl;
    // Check timing
    clock_t startTimer = clock();
    
    for (int k = 0; k < NTRY; k++) {
      MyCell *near = tree.getNearestNeighbour(xc[k]);
    }
    
    clock_t endTimer = clock();
    clock_t delta = endTimer-startTimer;
    double myTime = delta*1.0/CLOCKS_PER_SEC * 1000000.0;
    
    cout << "KDTree search/sec = " << myTime/NTRY << " us" << endl;
  
  }
  
  return 0;
}
