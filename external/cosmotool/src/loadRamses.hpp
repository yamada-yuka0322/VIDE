#ifndef _LOAD_RAMSES_HPP
#define _LOAD_RAMSES_HPP

#include "load_data.hpp"
#include "loadSimu.hpp"

namespace CosmoTool {

  GadgetData *loadRamses(const char *name, bool quiet = false);
  PurePositionData *loadRamsesPosition(const char *fname, int id, bool quiet = false, bool dp = true);
  PhaseSpaceData *loadRamsesPhase(const char *fname, int id, bool quiet = false);

  PhaseSpaceDataID *loadRamsesPhase1(const char *fname, int id, int cpuid, bool dp = true, bool quiet = false);

  SimuData *loadRamsesSimu(const char *basename, int id, int cpuid, bool dp, int flags);
};

#endif
