/*+
    VIDE -- Void IDentification and Examination -- ./c_tools/zobov2/voz1b1/voz_io.cpp
    Copyright (C) 2010-2014 Guilhem Lavaux
    Copyright (C) 2011-2014 P. M. Sutter

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; version 2 of the License.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
+*/
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
