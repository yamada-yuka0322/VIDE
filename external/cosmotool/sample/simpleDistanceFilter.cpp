#include "openmp.hpp"
#include "omptl/algorithm"
#include <cassert>
#include "yorick.hpp"
#include "sphSmooth.hpp"
#include "mykdtree.hpp"
#include "miniargs.hpp"
#include <H5Cpp.h>
#include "hdf5_array.hpp"
#include <iostream>
#include <boost/format.hpp>
#include <boost/bind.hpp>

using namespace std;
using namespace CosmoTool;

struct VCoord{
};

using boost::format;
using boost::str;
typedef boost::multi_array<float, 2> array_type_2d;
typedef boost::multi_array<float, 1> array_type_1d;

typedef KDTree<3,float> MyTree;
typedef MyTree::Cell MyCell;

int main(int argc, char **argv)
{

  char *fname1, *fname2;

  MiniArgDesc args[] = {
    { "INPUT DATA1", &fname1, MINIARG_STRING },
    { 0, 0, MINIARG_NULL }
  };

  if (!parseMiniArgs(argc, argv, args))
    return 1;

  H5::H5File in_f(fname1, 0);
  H5::H5File out_f("distances.h5", H5F_ACC_TRUNC);
  array_type_2d v1_data;
  array_type_1d dist_data;
  uint32_t N1_points;
  
  hdf5_read_array(in_f, "particles", v1_data);
  assert(v1_data.shape()[1] == 7);

  N1_points = v1_data.shape()[0];
  
  cout << "Got " << N1_points << " in the first file." << endl;

  MyCell *allCells_1 = new MyCell[N1_points];
  
  cout << "Shuffling data in cells..." << endl;
#pragma omp parallel for schedule(static)
  for (int i = 0 ; i < N1_points; i++)
    {
      for (int j = 0; j < 3; j++)
        allCells_1[i].coord[j] = v1_data[i][j];
      allCells_1[i].active = true;
    }
  dist_data.resize(boost::extents[N1_points]);
  
  cout << "Building trees..." << endl;
  MyTree tree1(allCells_1, N1_points);
#pragma omp parallel
  {
    MyCell **foundCells = new MyCell *[2];

    #pragma omp for
    for (size_t i = 0; i < N1_points; i++) {
      double dists[2];

      tree1.getNearestNeighbours(allCells_1[i].coord, 2, foundCells, dists);
      dist_data[i] = dists[1];
    }

    delete[] foundCells;
  }
  
  hdf5_write_array(out_f, "distances", dist_data);
  
  return 0;
};
