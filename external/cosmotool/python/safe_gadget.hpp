#include "config.hpp"
#include "loadGadget.hpp"
#include <string>

static inline
CosmoTool::SimuData *loadGadgetMulti_safe(const std::string& fname, int flags, int gadgetFormat)
{
  try
    {
      return CosmoTool::loadGadgetMulti(fname.c_str(), -1, flags, gadgetFormat);
    }
  catch (const CosmoTool::Exception& e)
    {
      return 0;
    }
}


static inline
CosmoTool::SimuData **alloc_simudata(int n)
{
  return new CosmoTool::SimuData *[n];
}

static inline
void del_simudata(CosmoTool::SimuData **s)
{
  delete[] s;
}
