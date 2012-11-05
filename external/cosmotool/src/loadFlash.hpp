#ifndef __COSMO_LOAD_FLASH_HPP
#define __COSMO_LOAD_FLASH_HPP

#include "load_data.hpp"
#include "loadSimu.hpp"

namespace CosmoTool {
  
  SimuData *loadFlashMulti(const char *fname, int id, int flags);
  
};

#endif
