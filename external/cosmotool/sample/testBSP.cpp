#if 0
#include "bsp_simple.hpp"

int main(int argc, char **argv)
{
  CosmoTool::simple_bsp::BSP<int, double, 2> bsp;
  double p[5][2] = { { 0, 0}, {1, 0}, {1, 1}, {0, 1}, {0, 0} };

  for (int q = 0; q < 4; q++)
    bsp.insert(p+q, p+q+2, q);

  return 0;
}
#endif
int main() {}
