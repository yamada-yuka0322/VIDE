#include <limits>
#include <iostream>
#include <fstream>
#include "voz_io.hpp"
#include "CosmoTool/fortran.hpp"
#include <cmath>

using namespace CosmoTool;
using namespace std;

bool PositionData::readFrom(const string& fname)
{
  try
    {
      UnformattedRead f(fname.c_str());

      f.beginCheckpoint();
      np = f.readInt32();
      f.endCheckpoint();

      for (int j = 0; j < 3; j++)
        {
          xyz[j] = new float[np];

          f.beginCheckpoint();
          for (int p = 0; p < np; p++)
            xyz[j][p] = f.readReal32();
          f.endCheckpoint();
        }
    }
  catch (const NoSuchFileException& e)
    {
      return false;
    }
  catch (const InvalidUnformattedAccess& e)
    {
      return false;
    }
  catch (const EndOfFileException& e)
    {
      return false;
    }

  return true;
}

void PositionData::findExtrema()
{
  const float BF = std::numeric_limits<float>::max();

  for (int i = 0; i < 3; i++)
    {
      xyz_min[i] = BF;
      xyz_max[i] = -BF;
    }

  V0=1;
  for (int i = 0; i < 3; i++)
    {
      for (pid_t p = 0; p < np; p++)
        {
          xyz_min[i] = min(xyz_min[i], xyz[i][p]);
          xyz_max[i] = max(xyz_max[i], xyz[i][p]);
        }
      V0 *= (xyz_max[i]-xyz_min[i]);
    }
}
