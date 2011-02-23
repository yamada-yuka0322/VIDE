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
