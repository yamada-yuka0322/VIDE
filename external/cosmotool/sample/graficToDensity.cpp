#include <cmath>
#include <iostream>
#include <cstdlib>
#include <boost/multi_array.hpp>
#include <H5Cpp.h>
#include "hdf5_array.hpp"
#include "miniargs.hpp"
#include "fortran.hpp"

using namespace std;
using namespace CosmoTool;

//#define GRAFIC_GUILHEM

int main(int argc, char **argv)
{
  uint32_t res;
  char *fname;
  int id;

  MiniArgDesc desc[] = {
    { "GRAFIC", &fname, MINIARG_STRING },
    { 0, 0, MINIARG_NULL }
  };

  if (!parseMiniArgs(argc, argv, desc))
    return 1;

  UnformattedRead ur(fname);

  ur.beginCheckpoint();
  int32_t nx = ur.readInt32();
  int32_t ny = ur.readInt32();
  int32_t nz = ur.readInt32();
  float dx = ur.readReal32();
  float xo = ur.readReal32();
  float yo = ur.readReal32();
  float zo = ur.readReal32();
  float astart = ur.readReal32();
  float omega_m = ur.readReal32();
  float omega_nu = ur.readReal32();
  float h0 = ur.readReal32();
#ifdef GRAFIC_GUILHEM
  float w0 = ur.readReal32();
#endif
  ur.endCheckpoint();

  cout << "Grafic file: Nx=" << nx << " Ny=" << ny << " Nz=" << nz << endl;
  cout << "a_start = " << astart << endl;
  cout << "z_start = " << 1/astart - 1 << endl;
  cout << "L = " << nx*dx << endl;

  boost::multi_array<float, 3> density(boost::extents[nx][ny][nz]);

  for (int32_t iz = 0; iz < nz; iz++)
    {
      ur.beginCheckpoint();
      for (int32_t iy = 0; iy < ny; iy++)
	{
	  for (int32_t ix = 0; ix < nx; ix++)
	    {
	      density[ix][iy][iz] = ur.readReal32();
	    }
	}
      ur.endCheckpoint();
    }

  H5::H5File f("density.h5", H5F_ACC_TRUNC);
  hdf5_write_array(f, "density", density); 

  return 0;
}
