/*+
This is CosmoTool (./sample/testkd3.cpp) -- Copyright (C) Guilhem Lavaux (2007-2014)

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

#include <ctime>
#include <cassert>
#include <cstdlib>
#include <iostream>
#include <fstream>
#define __KD_TREE_NUMNODES
#include "mykdtree.hpp"
#include "kdtree_splitters.hpp"
#include <limits>

#define NTRY 30
#define NGB 24
#define ND 3

using namespace std;
using namespace CosmoTool;

typedef KDTree<ND,char,ComputePrecision,KD_homogeneous_cell_splitter<ND, char> > MyTree;
//typedef KDTree<ND,char,ComputePrecision > MyTree;
typedef KDCell<ND,char> MyCell;

static double periodic(double a)
{
  while (a < -0.5)
    a+=1;
  while (a > 0.5)
    a -= 1;
  return a;
}

MyCell *findNearest(MyTree::coords& xc, MyCell *cells, uint32_t Ncells)
{
    MyCell *near2 = 0;
    double R2 = INFINITY;
    for (int i = 0; i < Ncells; i++)
      {
	double d2 = 0;
	for (int j = 0; j < ND; j++)
	  {
	    double delta = periodic(xc[j]-cells[i].coord[j]);
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
  uint32_t Ncells = 3000000;
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

  tree.setPeriodic(true, 1.0);

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
  MyCell **ngb = new MyCell *[NGB];
  double *distances = new double[NGB];
  
  ofstream fngb("nearest.txt");
  for (int k = 0; k < NTRY; k++) {
    cout << "Seed = " << xc[k][0] << " " << xc[k][1] << " " << xc[k][2] << endl;
    tree.getNearestNeighbours(xc[k], NGB, ngb, distances);
    int last = -1;

    for (uint32_t i = 0; i < NGB; i++)
      {
	if (ngb[i] == 0)
	  continue;

	last = i;

	double d2 = 0;
	for (int l = 0; l < 3; l++)
	  d2 += ({double delta = periodic(xc[k][l] - ngb[i]->coord[l]); delta*delta;});
	fngb << ngb[i]->coord[0] << " " << ngb[i]->coord[1] << " " << ngb[i]->coord[2] << " " << sqrt(d2) << " " << sqrt(distances[i]) << endl;
      }
    fngb << endl << endl;
    double farther_dist = distances[last];
    for (uint32_t i = 0; i < Ncells; i++)
      {
         bool found = false;
         // If the points is not in the list, it means it is farther than the farthest  point
	 for (int j =0; j < NGB; j++)
          {
             if (&cells[i] == ngb[j]) {
               found = true;
               break;
             }
          } 
	 double dist_to_seed = 0;
	 for (int l = 0; l < 3; l++)
	   {
	     double delta = periodic(xc[k][l]-cells[i].coord[l]);
	     dist_to_seed += delta*delta; 
	   }
	 double delta_epsilon = fabs(farther_dist/dist_to_seed-1);
	 if (!found) 
          { 
	    if (dist_to_seed <= farther_dist &&
		delta_epsilon > 1000*sqrt(numeric_limits<double>::epsilon()))
	      {
		cout << "SHOULD BE IN (d=" << dist_to_seed << "): Failed iteration " << k << " particle " << i << " (farthest=" << sqrt(farther_dist) << ")" << endl;
		cout << "delta_epsilon = " << delta_epsilon << endl;
		abort();
	      }
          }
	 else
	  {	    
	    if (dist_to_seed > farther_dist && 
		delta_epsilon > 1000*sqrt(numeric_limits<double>::epsilon())) 
	      {	
		cout << "SHOULD BE OUT (d=" << sqrt(dist_to_seed) << "): Failed iteration " << k << " particle " << i << " (farthest=" << sqrt(farther_dist) << ")" << endl;
		cout << "delta_epsilon = " << delta_epsilon << " epsilon = " << 100*sqrt(numeric_limits<double>::epsilon()) << endl;
		abort();
	      }
	  }
      }
  }

  return 0;
}
