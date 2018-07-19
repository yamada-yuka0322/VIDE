/*+
This is CosmoTool (./sample/testkd.cpp) -- Copyright (C) Guilhem Lavaux (2007-2014)

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
