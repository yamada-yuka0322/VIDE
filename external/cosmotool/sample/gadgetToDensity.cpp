/*+
This is CosmoTool (./sample/gadgetToDensity.cpp) -- Copyright (C) Guilhem Lavaux (2007-2014)

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
#include <iostream>
#include <cstdlib>
#include "cic.hpp"
#include "loadGadget.hpp"
#include "miniargs.hpp"
#include "yorick.hpp"

using namespace std;
using namespace CosmoTool;

int main(int argc, char **argv)
{
  uint32_t res;
  char *fname;
  int id;

  MiniArgDesc desc[] = {
    { "SNAPSHOT", &fname, MINIARG_STRING },
    { "ID", &id, MINIARG_INT },
    { "RESOLUTION", &res, MINIARG_INT },
    { 0, 0, MINIARG_NULL }
  };

  if (!parseMiniArgs(argc, argv, desc))
    return 1;

  SimuData *p = loadGadgetMulti(fname, 0, 0);
  double L0 = p->BoxSize;
  CICFilter filter(res, L0);

  delete p;

  try {
  for (int cpuid=0;;cpuid++) {
    p = loadGadgetMulti(fname, cpuid, NEED_POSITION);
    for (uint32_t i = 0; i < p->NumPart; i++)
    {
      CICParticles a;
      
      a.mass = 1.0;
      a.coords[0] = p->Pos[0][i]/1000;
      a.coords[1] = p->Pos[1][i]/1000;
      a.coords[2] = p->Pos[2][i]/1000;
      filter.putParticles(&a, 1);
    }
  
    delete p;
  }
  } catch (const NoSuchFileException& e) {}

  CICType *denField;
  uint32_t Ntot;
  filter.getDensityField(denField, Ntot);
  
  cout << "L0=" << L0 << endl;
  cout << "Saving density field" << endl;
  uint32_t dimList[] = { res, res, res};
  saveArray("densityField.nc", denField, dimList, 3);

  return 0;
}
