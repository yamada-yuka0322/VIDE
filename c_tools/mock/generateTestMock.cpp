/*+
    VIDE -- Void IDEntification pipeline -- ./c_tools/mock/generateTestMock.cpp
    Copyright (C) 2010-2013 Guilhem Lavaux
    Copyright (C) 2011-2013 Paul M. Sutter

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
+*/

#include <iostream>
#include <cstdlib>
#include <CosmoTool/fortran.hpp>

using namespace CosmoTool;
using namespace std;

#define LX 1.0
#define LY 1.0
#define LZ 1.0

#define NUMPART (16*16*16)

int main(int argc, char **argv)
{
  UnformattedWrite f("particles.bin");
  
  f.beginCheckpoint();
  f.writeInt32(NUMPART);
  f.endCheckpoint();

  cout << "Writing X components..." << endl;
  f.beginCheckpoint();
  for (uint32_t i = 0; i < NUMPART; i++)
    {
      f.writeReal32(drand48()*LX);
    }
  f.endCheckpoint();

  cout << "Writing Y components..." << endl;
  f.beginCheckpoint();
  for (uint32_t i = 0; i < NUMPART; i++)
    {
      f.writeReal32(drand48()*LY);
    }
  f.endCheckpoint();
  
  cout << "Writing Z components..." << endl;
  f.beginCheckpoint();
  for (uint32_t i = 0; i < NUMPART; i++)
    {
      f.writeReal32(drand48()*LZ);
    }
  f.endCheckpoint();

  return 0;
}
