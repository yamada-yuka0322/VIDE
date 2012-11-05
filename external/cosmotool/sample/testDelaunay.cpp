#include <fstream>
#include <iostream>
#include "dinterpolate.hpp"

#define NX 30
#define NY 30

using namespace std;
using namespace CosmoTool;

typedef DelaunayInterpolate<double,double,2> myTriangle;

int main()
{
  myTriangle::CoordType pos[] = { {0,0}, {1,0}, {0,1}, {1, 1}, {2, 0}, {2, 1} } ;
  double vals[] = { 0, 1, 1, 0, -0.5, 2.0 };
  uint32_t simplex[] = { 0, 1, 2, 3, 1, 2, 1, 3, 4, 4, 5, 3 };

  myTriangle t(&pos[0], &vals[0], &simplex[0], 6, 4);
  ofstream f("output.txt");

  for (uint32_t iy = 0; iy <= NY; iy++) {
    for (uint32_t ix = 0; ix <= NX; ix++) {
      myTriangle::CoordType inter = { ix *2.0/ NX, iy *1.0/NY };
      f << inter[1] << " " << inter[0] << " " << t.computeValue(inter) << endl;
    }
    f << endl;
  }
  return 0;
}
