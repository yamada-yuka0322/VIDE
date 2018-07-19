/*+
This is CosmoTool (./sample/testSmooth.cpp) -- Copyright (C) Guilhem Lavaux (2007-2014)

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

#include <iostream>
#include "sphSmooth.hpp"
#include "yorick.hpp"
#include "mykdtree.hpp"

using namespace std;
using namespace CosmoTool;

#define NX 1024
#define ND 2

typedef SPHSmooth<char,ND> MySmooth;
typedef MySmooth::SPHTree MyTree;
typedef MySmooth::SPHNode MyNode;
typedef MySmooth::SPHCell MyCell;

double unit_fun(const char& c)
{
  return 1.0;
}

int main()
{
  uint32_t Ncells = 10000;
  MyCell *cells = new MyCell[Ncells];
  
  for (int i = 0; i < Ncells; i++)
    {
      cells[i].active = true;
      for (int l = 0; l < ND; l++)
	cells[i].coord[l] = drand48();
      cells[i].val.weight = 0;
    }
  
  MyTree tree(cells, Ncells);
  MySmooth smooth(&tree, 16);

  for (uint32_t iy = 0; iy < NX; iy++)
    {
      cout << "iy=" << iy << endl;
      for (uint32_t ix = 0; ix < NX; ix++)
        {
           MyTree::coords c = { 1.0f*ix/NX, 1.0f*iy/NX };
           smooth.fetchNeighbours(c);
           smooth.addGridSite(c);
        }
    }

  
  uint32_t dims[] = { NX, NX };
  ProgressiveOutput<ComputePrecision> out = 
    ProgressiveOutput<ComputePrecision>::saveArrayProgressive("out.nc", dims, 2);
//#pragma omp parallel for schedule(static)
  for (uint32_t iy = 0; iy < NX; iy++)
    {
      MySmooth::SPHState state;
      cout << "iy=" << iy << endl;
      for (uint32_t ix = 0; ix < NX; ix++)
        {
           MyTree::coords c = { 1.0f*ix/NX, 1.0f*iy/NX };
           smooth.fetchNeighbours(c, &state);
           out.put(smooth.computeSmoothedValue(c, unit_fun, &state));

	  
        }
    }


  return 0;
}
