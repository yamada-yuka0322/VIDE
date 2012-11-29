#ifndef __COSMO_LOAD_GADGET_HPP
#define __COSMO_LOAD_GADGET_HPP

#include "load_data.hpp"
#include "loadSimu.hpp"

namespace CosmoTool {

  SimuData *loadGadgetMulti(const char *fname, int id, int flags, 
			    int GadgetFormat = 1, SimuFilter filter = 0);
 
  // Only single snapshot supported
  void writeGadget(const char *fname, SimuData *data, int GadgetFormat = 1);
 
};

#endif
