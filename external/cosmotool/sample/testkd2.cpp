#include <ctime>
#include <cassert>
#include <cstdlib>
#include <iostream>
#include <fstream>
#define __KD_TREE_NUMNODES
#include "mykdtree.hpp"
#include "kdtree_splitters.hpp"

#define NTRY 10
#define ND 3

using namespace std;
using namespace CosmoTool;

typedef KDTree<ND,char,ComputePrecision,KD_homogeneous_cell_splitter<ND, char> > MyTree;
//typedef KDTree<ND,char,ComputePrecision > MyTree;
typedef KDCell<ND,char> MyCell;

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
	    R2 = d2;
	  }
      }
    return near2;
}

int main()
{
  uint32_t Ncells = 10000000;
  MyCell *cells = new MyCell[Ncells];
  
  for (int i = 0; i < Ncells; i++)
    {
      cells[i].active = true;
      for (int l = 0; l < ND; l++)
	cells[i].coord[l] = drand48();
    }
  
  // Check timing
  clock_t startTimer = clock();
  MyTree tree(cells, Ncells);
  clock_t endTimer = clock();

  clock_t delta = endTimer-startTimer;
  double myTime = delta*1.0/CLOCKS_PER_SEC * 1.0;

  cout << "KDTree build = " << myTime << " s" << endl;

  MyTree::coords *xc = new MyTree::coords[NTRY];

  cout << "Generating seeds..." << endl;
  for (int k = 0; k < NTRY; k++)
    {
      for (int l = 0; l < ND; l++)
	xc[k][l] = drand48();
    }

  // Check consistency
  cout << "Check consistency..." << endl;
  MyCell **ngb = new MyCell *[12];
  double *distances = new double[12];
  
  ofstream fngb("nearest.txt");
  for (int k = 0; k < NTRY; k++) {
    cout << "Seed = " << xc[k][0] << " " << xc[k][1] << " " << xc[k][2] << endl;
    tree.getNearestNeighbours(xc[k], 12, ngb, distances);
    int last = -1;

    for (uint32_t i = 0; i < 12; i++)
      {
	if (ngb[i] == 0)
	  continue;

	last = i;

	double d2 = 0;
	for (int l = 0; l < 3; l++)
	  d2 += ({double delta = xc[k][l] - ngb[i]->coord[l]; delta*delta;});
	fngb << ngb[i]->coord[0] << " " << ngb[i]->coord[1] << " " << ngb[i]->coord[2] << " " << sqrt(d2) << endl;
      }
    fngb << endl << endl;
    double farther_dist = distances[last];
    for (uint32_t i = 0; i < Ncells; i++)
      {
         bool found = false;
         // If the points is not in the list, it means it is farther than the farthest  point
	 for (int j =0; j < 12; j++)
          {
             if (&cells[i] == ngb[j]) {
               found = true;
               break;
             }
          } 
	 double dist_to_seed = 0;
	 for (int l = 0; l < 3; l++)
	   {
	     double delta = xc[k][l]-cells[i].coord[l];
	     dist_to_seed += delta*delta; 
	   }
        if (!found) 
          { 
	    if (dist_to_seed <= farther_dist)
	      abort();
          }
	else
	  {
	    if (dist_to_seed > farther_dist)
		abort();
	  }
      }
  }

  return 0;
}
