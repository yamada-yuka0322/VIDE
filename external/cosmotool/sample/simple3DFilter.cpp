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

#define N_SPH 32

struct VCoord{
  float v[3];
  float mass;
};

using boost::format;
using boost::str;
typedef boost::multi_array<float, 2> array_type;
typedef boost::multi_array<float, 3> array3_type;
typedef boost::multi_array<float, 4> array4_type;

ComputePrecision getVelocity(const VCoord& v, int i)
{
  return v.mass * v.v[i];
}

ComputePrecision getMass(const VCoord& v)
{
  return v.mass;
}

typedef SPHSmooth<VCoord> MySmooth;
typedef MySmooth::SPHTree MyTree;
typedef MyTree::Cell MyCell;

template<typename FuncT>
void computeInterpolatedField(MyTree *tree1, double boxsize, int Nres, double cx, double cy, double cz,
                              array3_type& bins, array3_type& arr, FuncT func, double rLimit2)
{
#pragma omp parallel
  {
    MySmooth smooth1(tree1, N_SPH);
    
#pragma omp for schedule(dynamic) 
    for (int rz = 0; rz < Nres; rz++)
      {
        double pz = (rz)*boxsize/Nres-cz;

        cout << format("[%d] %d / %d") % smp_get_thread_id() % rz % Nres << endl;
        for (int ry = 0; ry < Nres; ry++)
          {
            double py = (ry)*boxsize/Nres-cy;
            for (int rx = 0; rx < Nres; rx++)
              {
                double px = (rx)*boxsize/Nres-cx;
                
                MyTree::coords c = { float(px), float(py), float(pz) };

                double r2 = c[0]*c[0]+c[1]*c[1]+c[2]*c[2];
                if (r2 > rLimit2)
                  {
                    arr[rx][ry][rz] = 0;
                    continue;
                  }

                uint32_t numInCell = bins[rx][ry][rz];
                if (numInCell > N_SPH)
                  smooth1.fetchNeighbours(c, numInCell);
                else
                  smooth1.fetchNeighbours(c);

                arr[rx][ry][rz] = smooth1.computeSmoothedValue(c, func);
              }
          }
      }
  }
}

int main(int argc, char **argv)
{

  char *fname1, *fname2;
  double rLimit, boxsize, rLimit2, cx, cy, cz;
  int Nres;

  MiniArgDesc args[] = {
    { "INPUT DATA1", &fname1, MINIARG_STRING },
    { "RADIUS LIMIT", &rLimit, MINIARG_DOUBLE },
    { "BOXSIZE", &boxsize, MINIARG_DOUBLE },
    { "RESOLUTION", &Nres, MINIARG_INT },
    { "CX", &cx, MINIARG_DOUBLE },
    { "CY", &cy, MINIARG_DOUBLE },
    { "CZ", &cz, MINIARG_DOUBLE },
    { 0, 0, MINIARG_NULL }
  };

  if (!parseMiniArgs(argc, argv, args))
    return 1;

  H5::H5File in_f(fname1, 0);
  H5::H5File out_f("fields.h5", H5F_ACC_TRUNC);
  array_type v1_data;
  uint32_t N1_points, N2_points;
  
  array3_type bins(boost::extents[Nres][Nres][Nres]);

  rLimit2 = rLimit*rLimit;

  hdf5_read_array(in_f, "particles", v1_data);
  assert(v1_data.shape()[1] == 7);

  N1_points = v1_data.shape()[0];
  
  cout << "Got " << N1_points << " in the first file." << endl;

  MyCell *allCells_1 = new MyCell[N1_points];
  
#pragma omp parallel for schedule(static)
  for (long i = 0; i < Nres*Nres*Nres; i++)
    bins.data()[i] = 0;

  cout << "Shuffling data in cells..." << endl;
#pragma omp parallel for schedule(static)
  for (int i = 0 ; i < N1_points; i++)
    {
      for (int j = 0; j < 3; j++)
        allCells_1[i].coord[j] = v1_data[i][j];
      for (int k = 0; k < 3; k++)
        allCells_1[i].val.pValue.v[k] = v1_data[i][3+k];
      allCells_1[i].val.pValue.mass = v1_data[i][6];
      allCells_1[i].active = true;
      allCells_1[i].val.weight = 0.0;

      long rx = floor((allCells_1[i].coord[0]+cx)*Nres/boxsize+0.5);
      long ry = floor((allCells_1[i].coord[1]+cy)*Nres/boxsize+0.5);
      long rz = floor((allCells_1[i].coord[2]+cz)*Nres/boxsize+0.5);
      
      if (rx < 0 || rx >= Nres || ry < 0 || ry >= Nres || rz < 0 || rz >= Nres)
        continue;
	  
#pragma omp atomic update
      bins[rx][ry][rz]++;
    }
  v1_data.resize(boost::extents[1][1]);
  
  hdf5_write_array(out_f, "num_in_cell", bins);

  cout << "Building trees..." << endl;
  MyTree tree1(allCells_1, N1_points);

  cout << "Creating smoothing filter..." << endl;

//  array3_type out_rad_1(boost::extents[Nres][Nres][Nres]);
 
  cout << "Weighing..." << endl;

#pragma omp parallel
  {
    MySmooth smooth1(&tree1, N_SPH);
    
#pragma omp for schedule(dynamic) 
    for (int rz = 0; rz < Nres; rz++)
     {
        double pz = (rz)*boxsize/Nres-cz;

        (cout << rz << " / " << Nres << endl).flush();
        for (int ry = 0; ry < Nres; ry++)
          {
            double py = (ry)*boxsize/Nres-cy;
            for (int rx = 0; rx < Nres; rx++)
              {
                double px = (rx)*boxsize/Nres-cx;
                
                MyTree::coords c = { float(px), float(py), float(pz) };

                double r2 = c[0]*c[0]+c[1]*c[1]+c[2]*c[2];
                if (r2 > rLimit2)
                  {
                    continue;
                  }

                uint32_t numInCell = bins[rx][ry][rz];
                if (numInCell > N_SPH)
                  smooth1.fetchNeighbours(c, numInCell);
                else
                  smooth1.fetchNeighbours(c);
#pragma omp critical
                smooth1.addGridSite(c);
              }
          }
        (cout << " Done " << rz << endl).flush();
      }  
  }
  
  cout << "Interpolating..." << endl;

  array3_type interpolated(boost::extents[Nres][Nres][Nres]);
  
  computeInterpolatedField(&tree1, boxsize, Nres, cx, cy, cz,
                           bins, interpolated, getMass, rLimit2);
  hdf5_write_array(out_f, "density", interpolated);
  //out_f.flush();
  for (int i = 0; i < 3; i++) {
      computeInterpolatedField(&tree1, boxsize, Nres, cx, cy, cz,
                               bins, interpolated, boost::bind(getVelocity, _1, i), rLimit2);
      hdf5_write_array(out_f, str(format("p%d") % i), interpolated);
  }
  
  return 0;
};
