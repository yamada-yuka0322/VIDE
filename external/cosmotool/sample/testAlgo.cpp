#include <iostream>
#include <stdio.h>
#include "algo.hpp"

using namespace CosmoTool;
using namespace std;

int main(int argc, char **argv)
{
  cout << square(2) << endl;
  cout << cube(2) << endl;
  cout << square(2.1f) << endl;
  cout << cube(2.1f) << endl;
  
  cout << spower<-2>(2.1f) << endl;
  cout << spower<2>(2.1f) << endl;
  cout << spower<3>(2.1f) << endl;
  cout << spower<4>(2.1f) << endl;
  cout << spower<5>(2.1f) << endl;
  return 0;
}
