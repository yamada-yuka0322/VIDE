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
           MyTree::coords c = { 1.0*ix/NX, 1.0*iy/NX };
           smooth.fetchNeighbours(c);
           smooth.addGridSite(c);
        }
    }

  
  uint32_t dims[] = { NX, NX };
  ProgressiveOutput<ComputePrecision> out = 
    ProgressiveOutput<ComputePrecision>::saveArrayProgressive("out.nc", dims, 2);
  for (uint32_t iy = 0; iy < NX; iy++)
    {
      cout << "iy=" << iy << endl;
      for (uint32_t ix = 0; ix < NX; ix++)
        {
           MyTree::coords c = { 1.0*ix/NX, 1.0*iy/NX };
           smooth.fetchNeighbours(c);
           out.put(smooth.computeSmoothedValue(c, unit_fun));

	  
        }
    }


  return 0;
}
