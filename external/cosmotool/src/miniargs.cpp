/*+
This is CosmoTool (./src/miniargs.cpp) -- Copyright (C) Guilhem Lavaux (2007-2014)

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

#include <cstdlib>
#include <cstdio>
#include <cstring>
#include "miniargs.hpp"
#include <iostream>

using namespace CosmoTool;
using namespace std;

int CosmoTool::parseMiniArgs(int argc, char **argv, MiniArgDesc *desc)
{
  int numMini;
  for (numMini = 0; desc[numMini].name != 0; numMini++);

  if ((argc-1) != numMini)
    {
      cerr << "Usage: ";
      for (int i = 0; i < numMini; i++)
	{
	  cerr << '[' << desc[i].name << "] ";
	}
      cerr << endl;
      return 0;
    }

  for (int i = 0; i < numMini; i++)
    {
      switch (desc[i].argType)
	{
	case MINIARG_STRING:
	  *((char **)desc[i].data) = strdup(argv[i+1]);
	  break;
	case MINIARG_INT:
	  *((int *)desc[i].data) = strtol(argv[i+1], NULL, 0);
	  break;
	case MINIARG_DOUBLE:
	  *((double *)desc[i].data) = strtod(argv[i+1], NULL);
	  break;
        case MINIARG_FLOAT:
	  *((float *)desc[i].data) = strtod(argv[i+1], NULL);
          break;
	case MINIARG_DOUBLE_3D_VECTOR:
	  {
	    double *d_array = (double *)(desc[i].data);

	    if (sscanf(argv[i+1], "(%lg,%lg,%lg)", &d_array[0], &d_array[1], &d_array[2]) != 3)
	      return 0;
	    break;
	  }
	}
    }

  return 1;
}

