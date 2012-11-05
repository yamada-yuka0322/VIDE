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

