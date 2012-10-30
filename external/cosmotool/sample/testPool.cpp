#include "pool.hpp"

using namespace CosmoTool;

int main(int argc, char **argv)
{
  MemoryPool<int> pool(1024);
  int **j = new int *[3000];

  for (int i = 0; i < 3000; i++)
    {
      j[i] = pool.alloc();
      j[i][0] = i;
    }

  pool.free_all();

  return 0;
}
