/*+
This is CosmoTool (./sample/gadgetToArray.cpp) -- Copyright (C) Guilhem Lavaux (2007-2014)

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
#include <cmath>
#include <algorithm>
#include <iostream>
#include <cstdlib>
#include "cic.hpp"
#include "loadGadget.hpp"
#include "miniargs.hpp"
#include <H5Cpp.h>
#include "hdf5_array.hpp"

using namespace CosmoTool;
using namespace std;

int main(int argc, char **argv)
{
  typedef boost::multi_array<float, 2> array_type;
  uint32_t res;
  char *fname;
  int id;
  double MPC;

  MiniArgDesc desc[] = {
    { "SNAPSHOT", &fname, MINIARG_STRING },
    { "MPC", &MPC, MINIARG_DOUBLE },
    { 0, 0, MINIARG_NULL }
  };

  if (!parseMiniArgs(argc, argv, desc))
    return 1;

  H5::H5File f("density.h5", H5F_ACC_TRUNC);
  

  SimuData *p = loadGadgetMulti(fname, 0, 0);
  double L0 = p->BoxSize/MPC;
  cout << "Will read " << p->TotalNumPart << " particles" << endl;
  array_type parts(boost::extents[p->TotalNumPart][7]);
  uint64_t q = 0;
  
  try {
    for (int cpuid=0;;cpuid++) {
      cout << "  = CPU " << cpuid << " = " << endl;
      p = loadGadgetMulti(fname, cpuid, NEED_POSITION|NEED_VELOCITY|NEED_MASS);
      cout << "  = DONE LOAD, COPYING IN PLACE" << endl;
      for (uint32_t i = 0; i < p->NumPart; i++)
        {
          for (int j = 0; j < 3; j++)
            {
              parts[q][j] = p->Pos[j][i]/MPC;
              while (parts[q][j] < 0) parts[q][j] += L0;
              while (parts[q][j] >= L0) parts[q][j] -= L0;
              parts[q][j] -= L0/2;
            }
          parts[q][3] = p->Vel[0][i];
          parts[q][4] = p->Vel[1][i];
          parts[q][5] = p->Vel[2][i];
          parts[q][6] = p->Mass[i];
          q++;
        }
      cout << "  = DONE (q=" << q << ")" << endl;

      delete p;
    }
  } catch (const NoSuchFileException& e) {}

  cout << " ++ WRITING ++" << endl;
  hdf5_write_array(f, "particles", parts);

  return 0;
}
